"""Disk handoff for per-fragment Casida results (MPI-safe, low peak memory on rank 0)."""
from __future__ import annotations

import os
import shutil
from typing import Any, Dict, List, Optional

import numpy as np


def fragment_casida_path(scratch_dir: str, frag_idx: int) -> str:
    return os.path.join(scratch_dir, f"casida_frag_{frag_idx:04d}.npz")


def _shared_scratch_root(base: Optional[str] = None) -> str:
    """Filesystem visible to all MPI ranks (not node-local ``/tmp``)."""
    if base:
        return base
    for key in ("EDFTPY_CASIDA_SCRATCH", "SLURM_SUBMIT_DIR"):
        val = os.environ.get(key)
        if val:
            return val
    return os.getcwd()


def get_casida_scratch_dir(comm, is_mpi: bool, base: Optional[str] = None) -> str:
    """Job-unique directory on a shared filesystem (all ranks use the same path)."""
    root = _shared_scratch_root(base)
    job = os.environ.get("SLURM_JOB_ID") or f"local_{os.getpid()}"
    path = os.path.join(root, ".edftpy_casida_scratch", str(job))
    if is_mpi:
        if comm.rank == 0:
            os.makedirs(path, exist_ok=True)
        comm.Barrier()
    else:
        os.makedirs(path, exist_ok=True)
    return path


def _n_transition_densities(rho) -> int:
    """Number of transitions whether ``rho`` is a list or stacked ndarray."""
    if rho is None:
        return 0
    if isinstance(rho, np.ndarray):
        if rho.size == 0:
            return 0
        if rho.ndim >= 4:
            return int(rho.shape[0])
        if rho.ndim == 3:
            return 1
        if rho.dtype == object:
            return len(rho)
    return len(rho)


def _rho_transition_to_stack(rho) -> tuple[np.ndarray, int]:
    """Normalize list or stacked-array ``rho_transition`` to ``(stack, n_trans)``."""
    if rho is None:
        return np.zeros((0,), dtype=float), 0
    if isinstance(rho, np.ndarray):
        if rho.size == 0:
            return np.zeros((0,), dtype=float), 0
        if rho.ndim >= 4:
            stack = np.asarray(rho, dtype=float)
            return stack, int(stack.shape[0])
        if rho.ndim == 3:
            stack = np.asarray(rho, dtype=float)[np.newaxis, ...]
            return stack, 1
        if rho.dtype == object:
            rho = [rho[i] for i in range(len(rho))]
        else:
            raise ValueError(f"unsupported rho_transition ndarray shape {rho.shape}")
    stack = np.stack([np.asarray(x) for x in rho], axis=0)
    return stack, int(stack.shape[0])


def casida_results_without_rho(results: Dict[str, Any]) -> Dict[str, Any]:
    """Lightweight copy for MPI gather (no transition-density grids)."""
    meta = {k: v for k, v in results.items() if k != "rho_transition"}
    meta["n_trans"] = _n_transition_densities(results.get("rho_transition"))
    return meta


def write_fragment_casida(path: str, results: Dict[str, Any]) -> None:
    """Write one fragment's Casida payload to a compressed ``.npz``."""
    rho_stack, n_trans = _rho_transition_to_stack(results.get("rho_transition"))

    payload = {
        "omega": np.asarray(results["omega"]),
        "Z": np.asarray(results.get("Z", results.get("eigenvectors"))),
        "n_trans": np.array(n_trans, dtype=np.int64),
        "rho_transition": rho_stack,
    }
    f = results.get("f", results.get("os_strength"))
    if f is not None:
        payload["f"] = np.asarray(f)
    dip = results.get("dip_tran")
    if dip is not None:
        payload["dip_tran"] = np.asarray(dip)

    tmp = path + ".tmp.npz"
    np.savez_compressed(tmp, **payload)
    os.replace(tmp, path)


def load_fragment_casida_meta(path: str) -> Dict[str, Any]:
    """Load energies, eigenvectors, dipoles — not ``rho_transition``."""
    with np.load(path, allow_pickle=False) as z:
        meta: Dict[str, Any] = {
            "omega": np.asarray(z["omega"]),
            "Z": np.asarray(z["Z"]),
            "eigenvectors": np.asarray(z["Z"]),
            "n_trans": int(np.asarray(z["n_trans"]).item()),
        }
        if "f" in z:
            meta["f"] = np.asarray(z["f"])
            meta["os_strength"] = meta["f"]
        if "dip_tran" in z:
            meta["dip_tran"] = np.asarray(z["dip_tran"])
    return meta


def load_fragment_rho_transition(path: str) -> List[np.ndarray]:
    """Load transition densities for one fragment (list of 3D arrays)."""
    with np.load(path, allow_pickle=False) as z:
        stack = np.asarray(z["rho_transition"])
    if stack.size == 0:
        return []
    return [stack[i] for i in range(stack.shape[0])]


def write_local_fragment_files(
    drivers,
    scratch_dir: str,
    comm,
    is_mpi: bool,
) -> List[Optional[str]]:
    """Each rank writes its subsystem Casida results; returns per-index paths (rank 0 merged)."""
    nsub = len(drivers)
    local = []

    for idx, driver in enumerate(drivers):
        if driver is None:
            continue
        res = getattr(driver, "casida_results", None)
        if res is None:
            continue
        path = fragment_casida_path(scratch_dir, idx)
        write_fragment_casida(path, res)
        driver.casida_results = casida_results_without_rho(res)
        local.append((idx, path))

    if not is_mpi:
        paths = [None] * nsub
        for idx, path in local:
            paths[idx] = path
        return paths

    gathered = comm.gather(local, root=0)
    if comm.rank != 0:
        return None

    paths = [None] * nsub
    for contributions in gathered:
        for idx, path in contributions:
            paths[idx] = path
    return paths


def _cross_fragment_matched_state_indices(
    fragment_results: List[Optional[Dict[str, Any]]],
    energy_thresh: float,
) -> List[Optional[List[int]]]:
    """State indices per fragment within ``energy_thresh`` (Hartree) of some other fragment."""
    n = len(fragment_results)
    matched: List[set] = [set() for _ in range(n)]
    for I in range(n):
        res_i = fragment_results[I]
        if res_i is None:
            continue
        omega_i = np.asarray(res_i["omega"], dtype=float)
        for i in range(len(omega_i)):
            for J in range(n):
                if J == I:
                    continue
                res_j = fragment_results[J]
                if res_j is None:
                    continue
                omega_j = np.asarray(res_j["omega"], dtype=float)
                if np.any(np.abs(omega_i[i] - omega_j) < energy_thresh):
                    matched[I].add(i)
                    break
    out: List[Optional[List[int]]] = []
    for s in matched:
        out.append(sorted(s) if s else None)
    return out


def _collapse_transition_densities_to_state_basis(
    Z_cols: np.ndarray,
    rho_stack: np.ndarray,
) -> tuple[List[np.ndarray], np.ndarray]:
    """Build one transition density per kept excitation from its eigenvector column.

    ``Z_cols`` is ``(n_trans, n_states)``; column ``k`` is eigenvector ``k`` in the
    occ–unocc transition basis. Removed excitations drop their columns; each kept
    column yields one grid for Pavanello coupling (length ``n_states``, not ``n_trans``).
    """
    n_trans, n_states = Z_cols.shape
    phi_states: List[np.ndarray] = []
    for k in range(n_states):
        phi_k = np.zeros_like(rho_stack[0], dtype=float)
        for ia in range(n_trans):
            phi_k = phi_k + Z_cols[ia, k] * rho_stack[ia]
        phi_states.append(phi_k)
    return phi_states, np.eye(n_states, dtype=float)


def reduce_one_fragment_casida(
    results: Dict[str, Any],
    state_indices: List[int],
    z_row_eps: float = 0.0,
    rho_path: Optional[str] = None,
) -> Dict[str, Any]:
    """Subset excitations; drop their eigenvector columns and collapse ``rho`` to state basis."""
    state_ix = np.asarray(state_indices, dtype=int)
    omega = np.asarray(results["omega"], dtype=float)
    Z = np.asarray(results.get("Z", results.get("eigenvectors")), dtype=float)
    n_states_full = len(omega)
    n_trans = Z.shape[0]

    omega_k = omega[state_ix]
    Z_cols = Z[:, state_ix]
    n_kept = len(state_ix)

    reduced: Dict[str, Any] = {
        "omega": omega_k,
        "n_trans": n_kept,
    }

    f = results.get("f", results.get("os_strength"))
    if f is not None:
        f = np.asarray(f, dtype=float)
        reduced["f"] = f[state_ix]
        reduced["os_strength"] = reduced["f"]

    dip = results.get("dip_tran")
    if dip is not None:
        dip = np.asarray(dip, dtype=float)
        if dip.shape[0] == n_states_full:
            reduced["dip_tran"] = dip[state_ix]
        elif dip.shape[0] == n_trans:
            reduced["dip_tran"] = np.real(Z_cols.T @ dip)

    rho = results.get("rho_transition")
    if rho_path:
        rho = load_fragment_rho_transition(rho_path)
    if rho is not None:
        rho_stack, _ = _rho_transition_to_stack(rho)
        phi_states, Z_states = _collapse_transition_densities_to_state_basis(
            np.real(Z_cols), rho_stack,
        )
        reduced["rho_transition"] = phi_states
        reduced["Z"] = Z_states
        reduced["eigenvectors"] = Z_states
    else:
        reduced["Z"] = Z_cols
        reduced["eigenvectors"] = Z_cols

    return reduced


def reduce_active_space(
    fragment_results: List[Optional[Dict[str, Any]]],
    energy_thresh: float,
    z_row_eps: float = 0.0,
    stream_paths: Optional[List[Optional[str]]] = None,
) -> List[Optional[Dict[str, Any]]]:
    """Cross-fragment energy window on ``omega``; shrink ``Z`` and ``rho_transition`` accordingly."""
    state_sets = _cross_fragment_matched_state_indices(
        fragment_results, energy_thresh,
    )
    reduced: List[Optional[Dict[str, Any]]] = []
    for idx, res in enumerate(fragment_results):
        if res is None:
            reduced.append(None)
            continue
        state_ix = state_sets[idx]
        if not state_ix:
            state_ix = list(range(len(np.asarray(res["omega"]))))
        path = None
        if stream_paths and idx < len(stream_paths):
            path = stream_paths[idx]
        reduced.append(
            reduce_one_fragment_casida(
                res, state_ix, z_row_eps=z_row_eps, rho_path=path,
            ),
        )
    return reduced


def uncoupled_excluded_states(
    fragment_results: List[Optional[Dict[str, Any]]],
    energy_thresh: float,
    *,
    recompute_f: bool = True,
) -> List[Optional[List[Dict[str, Any]]]]:
    """Per-fragment local Casida data for states dropped by active-space reduction.

    Uses the same cross-fragment ``|Δω|`` rule as :func:`reduce_active_space`.
    When no state matches the threshold on a fragment, that fragment keeps all
    states for coupling and returns an empty exclusion list.
    """
    state_sets = _cross_fragment_matched_state_indices(
        fragment_results, energy_thresh,
    )
    excluded_by_frag: List[Optional[List[Dict[str, Any]]]] = []
    for idx, res in enumerate(fragment_results):
        if res is None:
            excluded_by_frag.append(None)
            continue
        omega = np.asarray(res["omega"], dtype=float)
        n_states = len(omega)
        kept = state_sets[idx]
        if not kept:
            excluded_ix: List[int] = []
        else:
            kept_set = set(kept)
            excluded_ix = [i for i in range(n_states) if i not in kept_set]
        f = res.get("f", res.get("os_strength"))
        entries: List[Dict[str, Any]] = []
        for i in excluded_ix:
            entry: Dict[str, Any] = {
                "state_index": int(i),
                "omega": float(omega[i]),
            }
            if recompute_f:
                entry["f"] = recompute_fragment_oscillator_strength(res, i)
            elif f is not None:
                entry["f"] = float(np.asarray(f, dtype=float)[i])
            entries.append(entry)
        excluded_by_frag.append(entries)
    return excluded_by_frag


_HARTREE_TO_EV = 27.211386246


def recompute_fragment_oscillator_strength(
    fragment_res: Dict[str, Any],
    state_index: int,
) -> float:
    """``f_n = (2/3) omega_n |d_n|^2`` with fragment ``omega``, ``Z``, ``dip_tran``.

    Returns the raw Casida oscillator strength (not sum-normalized). The merged
    spectrum applies a single ``sum(f)=1`` scaling in
    :func:`merge_coupled_and_uncoupled_spectrum`.
    """
    omega = np.asarray(fragment_res["omega"], dtype=float)
    w = float(omega[state_index])
    if w < 1e-30:
        return 0.0
    mu = fragment_res.get("dip_tran")
    if mu is None:
        f = fragment_res.get("f", fragment_res.get("os_strength"))
        if f is not None:
            return float(np.asarray(f, dtype=float)[state_index])
        return float("nan")
    mu = np.asarray(mu, dtype=float)
    Z = np.asarray(fragment_res.get("Z", fragment_res.get("eigenvectors")), dtype=float)
    n_states = len(omega)
    if mu.shape[0] == n_states:
        d = mu[state_index]
    else:
        d = np.real(mu.T @ Z[:, state_index])
    return (2.0 / 3.0) * w * float(np.dot(d, d))


def write_uncoupled_excluded_txt(
    path: str,
    excluded_by_frag: List[Optional[List[Dict[str, Any]]]],
) -> str:
    """Write uncoupled fragment-local ``omega`` (Ha, eV) and ``f`` to a text file."""
    path = os.path.abspath(path)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    n_total = 0
    with open(path, "w") as f:
        f.write(
            "# subsystem  state_index  omega_Ha  omega_eV  oscillator_strength\n",
        )
        for frag_idx, entries in enumerate(excluded_by_frag):
            if not entries:
                continue
            for entry in entries:
                omega_ha = float(entry["omega"])
                omega_ev = omega_ha * _HARTREE_TO_EV
                fval = entry.get("f", float("nan"))
                f.write(
                    f"{frag_idx:5d}  {int(entry['state_index']):5d}  "
                    f"{omega_ha:14.8f}  {omega_ev:14.8f}  {float(fval):14.8f}\n",
                )
                n_total += 1
    if n_total == 0:
        with open(path, "a") as f:
            f.write("# (no uncoupled excluded states)\n")
    return path


def merge_coupled_and_uncoupled_spectrum(
    coupled: Optional[Dict[str, Any]],
    excluded_by_frag: Optional[List[Optional[List[Dict[str, Any]]]]] = None,
    *,
    fragment_results_full: Optional[List[Optional[Dict[str, Any]]]] = None,
    sort_by_energy: bool = True,
    normalize_fosc: bool = True,
) -> Dict[str, Any]:
    """Build the final excitation spectrum: coupled block + fragment-local excluded states.

    Coupled states use ``coupled['omega']`` and ``coupled['f']`` (from
    :func:`casidapy.subsystem_coupling.coupled_oscillator_strengths`). Uncoupled
    states use fragment ``omega`` with ``f`` from ``recompute_fragment_oscillator_strength``
    when ``fragment_results_full`` is passed (same ``(2/3) omega |d|^2`` as coupling,
    without per-fragment normalization of stored Casida ``f``).

    Returns keys ``omega_all``, ``f_all``, ``n_coupled``, ``n_uncoupled``, and
    provenance arrays ``is_coupled``, ``fragment_index``, ``state_index``.

    If ``normalize_fosc`` is true (default), ``f_all`` is scaled so
    ``sum(f_all) == 1`` after coupled and uncoupled entries are merged.
    """
    omega_list: List[float] = []
    f_list: List[float] = []
    is_coupled: List[bool] = []
    fragment_index: List[int] = []
    state_index: List[int] = []

    if coupled is not None:
        om = coupled.get("omega")
        if om is not None:
            om = np.asarray(om, dtype=float).ravel()
            ff = coupled.get("f")
            if ff is not None:
                ff = np.asarray(ff, dtype=float).ravel()
                if len(ff) != len(om):
                    raise ValueError(
                        f"len(coupled f)={len(ff)} != len(coupled omega)={len(om)}",
                    )
            else:
                ff = np.full(len(om), np.nan, dtype=float)
            for w, fv in zip(om, ff):
                omega_list.append(float(w))
                f_list.append(float(fv))
                is_coupled.append(True)
                fragment_index.append(-1)
                state_index.append(-1)

    n_unc = 0
    if excluded_by_frag:
        for frag_idx, entries in enumerate(excluded_by_frag):
            if not entries:
                continue
            for entry in entries:
                i = int(entry["state_index"])
                omega_list.append(float(entry["omega"]))
                if fragment_results_full is not None and frag_idx < len(
                    fragment_results_full,
                ):
                    res_full = fragment_results_full[frag_idx]
                    if res_full is not None:
                        fv = recompute_fragment_oscillator_strength(res_full, i)
                    else:
                        fv = float(entry.get("f", float("nan")))
                else:
                    fv = float(entry.get("f", float("nan")))
                f_list.append(fv)
                is_coupled.append(False)
                fragment_index.append(int(frag_idx))
                state_index.append(int(entry.get("state_index", -1)))
                n_unc += 1

    n_cpl = int(sum(is_coupled))
    n_unc = len(is_coupled) - n_cpl
    omega_all = np.asarray(omega_list, dtype=float)
    f_all = np.asarray(f_list, dtype=float)
    is_coupled_arr = np.asarray(is_coupled, dtype=bool)
    frag_ix = np.asarray(fragment_index, dtype=np.int32)
    state_ix = np.asarray(state_index, dtype=np.int32)

    if sort_by_energy and omega_all.size:
        order = np.argsort(omega_all)
        omega_all = omega_all[order]
        f_all = f_all[order]
        is_coupled_arr = is_coupled_arr[order]
        frag_ix = frag_ix[order]
        state_ix = state_ix[order]

    if normalize_fosc and f_all.size:
        f_sum = float(np.nansum(f_all))
        if f_sum > 0.0 and np.isfinite(f_sum):
            f_all = f_all / f_sum

    return {
        "omega_all": omega_all,
        "f_all": f_all,
        "n_coupled": n_cpl,
        "n_uncoupled": n_unc,
        "is_coupled": is_coupled_arr,
        "fragment_index": frag_ix,
        "state_index": state_ix,
    }


def build_fragment_results_from_files(
    paths: List[Optional[str]],
) -> List[Optional[Dict[str, Any]]]:
    """Build ``fragment_results`` metadata on rank 0 (no ``rho_transition`` in memory)."""
    out: List[Optional[Dict[str, Any]]] = []
    for path in paths:
        if path is None or not os.path.isfile(path):
            out.append(None)
            continue
        out.append(load_fragment_casida_meta(path))
    return out


def cleanup_fragment_files(
    paths: List[Optional[str]],
    scratch_dir: Optional[str] = None,
    comm=None,
    is_mpi: bool = False,
    *,
    root_rank_only: bool = True,
) -> None:
    """Remove per-fragment ``.npz`` files and optionally the scratch directory.

    With MPI, only world rank 0 should delete (other ranks must not ``rmtree``
    while rank 0 is still reading files during coupling).
    """
    if is_mpi and comm is not None:
        comm.Barrier()

    if root_rank_only and is_mpi and comm is not None and comm.rank != 0:
        if is_mpi and comm is not None:
            comm.Barrier()
        return

    if paths:
        for path in paths:
            if path and os.path.isfile(path):
                try:
                    os.remove(path)
                except OSError:
                    pass

    if scratch_dir and os.path.isdir(scratch_dir):
        try:
            shutil.rmtree(scratch_dir, ignore_errors=True)
        except OSError:
            pass

    if is_mpi and comm is not None:
        comm.Barrier()
