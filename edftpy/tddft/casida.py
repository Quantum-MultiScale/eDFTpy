import os
import time
from contextlib import contextmanager
from typing import Optional

import numpy as np
from dftpy.mpi import MP, SerialComm

from edftpy.mpi import graphtopo, sprint
from edftpy.optimizer import Optimization
from edftpy.tddft.casida_fragment_io import (
    build_fragment_results_from_files,
    cleanup_fragment_files,
    get_casida_scratch_dir,
    merge_coupled_and_uncoupled_spectrum,
    reduce_active_space,
    uncoupled_excluded_states,
    write_local_fragment_files,
    write_uncoupled_excluded_txt,
)

_HARTREE_TO_EV = 27.211386
from edftpy.utils.common import Field, Grid


def _mpi_context(gsystem):
    """Return (graphtopo, MPI_COMM_WORLD communicator, is_mpi flag)."""
    gt = gsystem.graphtopo if gsystem is not None else None
    comm = gt.comm if gt is not None else graphtopo.comm
    is_mpi = bool(gt.is_mpi) if gt is not None else False
    return gt, comm, is_mpi


def _gather_subsystem_list(drivers, comm, is_mpi, nsub, getter):
    """Collect per-rank subsystem data into ``[None] * nsub`` on world rank 0.

    Each rank contributes ``(subsystem_index, value)`` pairs from drivers it owns.
    Uses ``gather`` (not ``allgather``) so large payloads such as ``casida_results``
    with full transition-density lists are not replicated on every rank.
    Non-root ranks return ``None``; coupling only needs the merged list on rank 0.
    """
    local = []
    for idx, driver in enumerate(drivers):
        if driver is None:
            continue
        value = getter(driver)
        if value is not None:
            local.append((idx, value))

    if not is_mpi:
        out = [None] * nsub
        for idx, value in local:
            out[idx] = value
        return out

    gathered = comm.gather(local, root=0)
    if comm.rank != 0:
        return None

    merged = [None] * nsub
    for contributions in gathered:
        if not contributions:
            continue
        for idx, value in contributions:
            merged[idx] = value
    return merged


def _prepare_serial_coupling_fields(gsystem, comm, is_mpi):
    """Build a serial global grid and total density for Pavanello coupling on rank 0.

    All ranks must call this: ``Field.gather`` is collective. Non-root ranks
    return ``(None, None)``.
    """
    rho = gsystem.density
    if getattr(rho, "rank", 1) > 1:
        rho = Field(
            grid=rho.grid,
            rank=1,
            data=np.sum(np.asarray(rho), axis=0),
        )

    rho_gathered = rho.gather(root=0) if is_mpi else rho
    if is_mpi and comm.rank != 0:
        return None, None

    mp_serial = MP(comm=SerialComm())
    g = gsystem.grid
    serial_grid = Grid(
        lattice=g.lattice,
        nr=g.nrR,
        full=g.full,
        mp=mp_serial,
        direct=True,
    )
    serial_rho = Field(
        grid=serial_grid,
        rank=1,
        data=np.asarray(rho_gathered),
    )
    return serial_grid, serial_rho


@contextmanager
def _serial_coupling_gsystem(gsystem, serial_grid, serial_rho):
    """Temporarily use gathered serial grid/density for Pavanello coupling on rank 0."""
    gt = gsystem.graphtopo
    saved = {
        "grid": gsystem._grid,
        "density": gsystem._density,
        "graph_grid": gt.graph.grid,
    }
    gsystem._grid = serial_grid
    gsystem._density = serial_rho
    gt.graph.grid = serial_grid
    try:
        yield gsystem
    finally:
        gsystem._grid = saved["grid"]
        gsystem._density = saved["density"]
        gt.graph.grid = saved["graph_grid"]


def _run_local_fragment_casida(drivers, options, comm, is_mpi):
    """Run Casida on every subsystem driver owned by this MPI rank.

    Fragment teams use ``comm_sub`` inside ``update_workspace`` (QEpy + CasidaPy).
    Different fragments proceed at the same wall-clock time; a single WORLD
    ``Barrier`` at the end synchronizes before global gather/coupling.
    """
    for idx, driver in enumerate(drivers):
        if driver is None:
            continue
        sprint(
            f"Fragment Casida subsystem {idx} ({driver.prefix})...",
            flush=True,
        )
        driver.save(save=["W", "D"])
        driver.task = "casida"
        driver.update_workspace(first=True, options=options)
        sprint(f"Fragment Casida subsystem {idx} done.", flush=True)

    if is_mpi:
        comm.Barrier()
        sprint("Fragment Casida: WORLD sync before gather/coupling.", flush=True)


def _slim_coupled_results(coupled):
    """MPI-broadcast payload: spectra only (omit large ``K_coupling`` / ``f_nadd``)."""
    if coupled is None:
        return None
    return {
        "omega": coupled.get("omega"),
        "f": coupled.get("f"),
        "Z": coupled.get("Z"),
        "omega_all": coupled.get("omega_all"),
        "f_all": coupled.get("f_all"),
        "n_coupled": coupled.get("n_coupled"),
        "n_uncoupled": coupled.get("n_uncoupled"),
        "is_coupled": coupled.get("is_coupled"),
        "fragment_index": coupled.get("fragment_index"),
        "state_index": coupled.get("state_index"),
        "plot_path": coupled.get("plot_path"),
    }


def _log_fragment_casida_results(fragment_results, olevel):
    """Print per-fragment Casida summaries; full arrays only if ``olevel >= 2``."""
    for idx, res in enumerate(fragment_results):
        if res is None:
            continue
        omega = res.get("omega")
        n_states = res.get("n_states", len(res.get("omega", [])))
        n_prim = res.get("n_trans_primitive", res.get("n_trans"))
        sprint(f"Subsystem {idx} Casida: {n_states} amplitude-basis rho grids "
               f"({n_prim} primitive transitions).")
        if olevel < 2:
            continue
        for key, value in res.items():
            if key == "rho_transition":
                rho_basis = res.get("rho_basis", "amplitude_xpy")
                sprint(
                    f"Subsystem {idx} Casida: {key}: {n_states} stored "
                    f"({rho_basis}; {n_prim} primitive transitions)",
                )
            else:
                sprint(f"Subsystem {idx} Casida: {key}:\n", value)


def _log_active_space_reduction(fragment_results, reduced_results, energy_thresh):
    """Summarize state/transition counts before and after reduction."""
    for idx, (before, after) in enumerate(zip(fragment_results, reduced_results)):
        if before is None or after is None:
            continue
        n_s0 = len(np.asarray(before["omega"]))
        n_s1 = len(np.asarray(after["omega"]))
        n_t0 = before.get("n_trans")
        if n_t0 is None:
            n_t0 = int(np.asarray(before.get("Z", [])).shape[0])
        n_t1 = after.get("n_trans", int(np.asarray(after.get("Z", [])).shape[0]))
        sprint(
            f"Subsystem {idx} active-space reduction (|Δω| < {energy_thresh} Ha): "
            f"{n_s0} -> {n_s1} states, {n_t0} -> {n_t1} transitions",
            flush=True,
        )


def _log_uncoupled_excluded_states(excluded_by_frag):
    """Print fragment-local energies/strengths for states not in the coupled active space."""
    n_total = sum(len(e) for e in excluded_by_frag if e)
    if n_total == 0:
        sprint(
            "No uncoupled excluded states (all fragment states entered coupling).",
            flush=True,
        )
        return
    sprint(
        f"Uncoupled excluded states ({n_total} total; fragment-local Casida values):",
        flush=True,
    )
    for frag_idx, entries in enumerate(excluded_by_frag):
        if not entries:
            continue
        sprint(f"  Subsystem {frag_idx}:", flush=True)
        for entry in entries:
            i = entry["state_index"]
            omega_ha = entry["omega"]
            omega_ev = omega_ha * _HARTREE_TO_EV
            line = (
                f"    state {i}: omega = {omega_ev:.6f} eV "
                f"({omega_ha:.6f} Ha)"
            )
            if "f" in entry:
                line += f", f = {entry['f']:.6f}"
            sprint(line, flush=True)


def _log_merged_spectrum(merged) -> None:
    """Print the combined coupled + uncoupled excitation spectrum (stdout + sprint)."""
    if not merged:
        return
    omega = np.asarray(merged.get("omega_all", []), dtype=float)
    f = np.asarray(merged.get("f_all", []), dtype=float)
    n_cpl = int(merged.get("n_coupled", 0))
    n_unc = int(merged.get("n_uncoupled", 0))
    header = (
        f"Merged Casida spectrum: {len(omega)} states "
        f"({n_cpl} coupled + {n_unc} uncoupled excluded), sorted by energy."
    )
    print(header, flush=True)
    sprint(header, flush=True)
    if omega.size == 0:
        return
    omega_ev = omega * _HARTREE_TO_EV
    print("Merged excitation energies (eV):", flush=True)
    print(omega_ev, flush=True)
    print("Merged oscillator strengths:", flush=True)
    print(f, flush=True)
    sprint(f"Merged excitation energies (eV):\n{omega_ev}", flush=True)
    sprint(f"Merged oscillator strengths:\n{f}", flush=True)
    is_cpl = merged.get("is_coupled")
    for i, (w, fv) in enumerate(zip(omega, f)):
        tag = "coupled"
        if is_cpl is not None and not np.asarray(is_cpl)[i]:
            tag = "uncoupled"
        line = (
            f"  [{i:4d}] {tag:9s}  omega = {w * _HARTREE_TO_EV:10.4f} eV "
            f"({w:.6f} Ha)"
        )
        if np.isfinite(fv):
            line += f",  f = {fv:.6f}"
        print(line, flush=True)
        sprint(line, flush=True)


def _plot_merged_spectrum(coupled, options) -> Optional[str]:
    """Write merged spectrum PNG via :mod:`casidapy.plot_casida_spectrum` (rank 0 only)."""
    if not options.get("casida_plot", True):
        return None
    omega = coupled.get("omega_all")
    f = coupled.get("f_all")
    if omega is None or f is None:
        return None
    omega = np.asarray(omega, dtype=float).ravel()
    f = np.asarray(f, dtype=float).ravel()
    if omega.size == 0:
        return None

    out_png = options.get("casida_plot_out") or "casida_merged_spectrum.png"
    if not os.path.isabs(out_png):
        root = os.environ.get("SLURM_SUBMIT_DIR") or os.getcwd()
        out_png = os.path.join(root, out_png)

    txt_name = options.get("casida_plot_txt") or "casida_merged_spectrum.txt"
    if not os.path.isabs(txt_name):
        root = os.environ.get("SLURM_SUBMIT_DIR") or os.getcwd()
        txt_name = os.path.join(root, txt_name)

    is_coupled = coupled.get("is_coupled")
    energies_ev = omega * _HARTREE_TO_EV

    try:
        from casidapy.plot_casida_spectrum import (
            plot_casida_spectrum,
            write_casida_spectrum_txt,
        )

        write_casida_spectrum_txt(txt_name, energies_ev, f, is_coupled=is_coupled)
        png_path = plot_casida_spectrum(
            energies_ev,
            f,
            out=out_png,
            sigma=float(options.get("casida_plot_sigma", 0.08)),
            is_coupled=is_coupled,
            title="Merged Casida spectrum (coupled + uncoupled)",
        )
        msg = f"Casida spectrum plot: {png_path}  (data: {txt_name})"
        print(msg, flush=True)
        sprint(msg, flush=True)
        return png_path
    except Exception as exc:
        warn = f"WARNING: Casida spectrum plot failed: {exc}"
        print(warn, flush=True)
        sprint(warn, flush=True)
        return None


class CasidaTDDFT(Optimization):
    """
    Subsystem Casida LR-TDDFT with inter-fragment coupling.

    MPI layout:
      - Embedded SCF: parallel on graphtopo.comm_sub per fragment.
      - Fragment Casida: each rank runs its own subsystem on comm_sub; fragments
        overlap in time (disjoint QE communicators). One WORLD Barrier before gather.
      - Coupling: transition densities file-streamed to rank 0 (``.npz`` per fragment).
    """

    def __init__(self, drivers=None, gsystem=None, options=None, optimizer=None, **kwargs):
        super().__init__(drivers=drivers, gsystem=gsystem, options=options)
        default_options = {
            "olevel": 1,
            "sdft": "sdft",
            "number_of_states": 10,
            "number_of_bands": None,
            "coupling": True,
            "tda": False,
            "rho_cutoff": 1e-3,
            "fxc_max": 20.0,
            "matrix_free": True,
            "solver_method": "eigsh",
            "coupling_energy_threshold": None,
            "coupling_z_row_eps": 0.0,
            "casida_plot": True,
            "casida_plot_sigma": 0.08,
            "casida_plot_out": "casida_merged_spectrum.png",
            "casida_plot_txt": "casida_merged_spectrum.txt",
            "casida_uncoupled_txt": "casida_uncoupled_spectrum.txt",
        }
        self.options = default_options
        if isinstance(options, dict):
            self.options.update(options)
        self.optimizer = optimizer
        self.coupled_results = None
        self._uncoupled_excluded = None

    def initialization(self):
        """Embedded SCF, then parallel per-fragment Casida on comm_sub (WORLD sync once)."""
        self.optimizer.optimize()
        _, comm, is_mpi = _mpi_context(self.gsystem)
        _run_local_fragment_casida(self.drivers, self.options, comm, is_mpi)

    def optimize(self, **kwargs):
        """Alias for :meth:`run` (Optimization API)."""
        self.run(**kwargs)

    def run(self, **kwargs):
        """Consume :meth:`irun` and return the last yielded result."""
        for item in self.irun(**kwargs):
            pass
        return item

    def irun(self, restart=None, **kwargs):
        """Run fragment Casida, gather results, optional Pavanello coupling, yield spectrum."""
        self.initialization()
        self.time_begin = time.time()

        _, comm, is_mpi = _mpi_context(self.gsystem)
        nsub = len(self.drivers)
        use_coupling = self.options.get("coupling", True)
        self._fragment_stream_paths = None
        self._fragment_disk_paths = None
        self._casida_scratch_dir = None
        self._uncoupled_excluded = None
        self._fragment_results_full = None
        self.gsystem.casida_uncoupled_excluded = None
        self.gsystem.casida_merged_spectrum = None

        if use_coupling:
            self._casida_scratch_dir = get_casida_scratch_dir(comm, is_mpi)
            if (not is_mpi) or comm.rank == 0:
                sprint(
                    f"Casida file stream directory: {self._casida_scratch_dir}",
                    flush=True,
                )
            self._fragment_stream_paths = write_local_fragment_files(
                self.drivers,
                self._casida_scratch_dir,
                comm,
                is_mpi,
            )
            if is_mpi:
                comm.Barrier()
            if (not is_mpi) or comm.rank == 0:
                paths = self._fragment_stream_paths or []
                missing = [
                    i for i, p in enumerate(paths) if p and not os.path.isfile(p)
                ]
                if missing:
                    sprint(
                        f"ERROR: Casida stream files missing on rank 0 for "
                        f"subsystems {missing} (use a shared directory, not "
                        f"node-local /tmp).",
                        flush=True,
                    )
                meta_before = build_fragment_results_from_files(paths)
                self._fragment_disk_paths = paths
                self._fragment_stream_paths = paths
                energy_thresh = self.options.get("coupling_energy_threshold")
                if energy_thresh is not None and float(energy_thresh) > 0:
                    z_row_eps = float(
                        self.options.get("coupling_z_row_eps", 0.0),
                    )
                    sprint(
                        f"Reducing Casida active space (|Δω| < "
                        f"{float(energy_thresh)} Ha)...",
                        flush=True,
                    )
                    eth = float(energy_thresh)
                    reduced = reduce_active_space(
                        meta_before,
                        energy_thresh=eth,
                        z_row_eps=z_row_eps,
                        stream_paths=paths,
                    )
                    _log_active_space_reduction(meta_before, reduced, eth)
                    self._fragment_results_full = meta_before
                    self._uncoupled_excluded = uncoupled_excluded_states(
                        meta_before, eth, recompute_f=True,
                    )
                    self.gsystem.casida_uncoupled_excluded = (
                        self._uncoupled_excluded
                    )
                    self.gsystem.casida_results = reduced
                    self._fragment_stream_paths = None
                else:
                    self.gsystem.casida_results = meta_before
            else:
                self.gsystem.casida_results = None
            if is_mpi:
                comm.Barrier()
        else:
            self.gsystem.casida_results = _gather_subsystem_list(
                self.drivers,
                comm,
                is_mpi,
                nsub,
                lambda d: getattr(d, "casida_results", None),
            )

        if (not is_mpi) or comm.rank == 0:
            n_ok = sum(1 for r in (self.gsystem.casida_results or []) if r is not None)
            sprint(
                f"Casida metadata for {n_ok}/{nsub} subsystems on rank 0.",
                flush=True,
            )
            _log_fragment_casida_results(
                self.gsystem.casida_results,
                self.options.get("olevel", 1),
            )

        sprint(f"Per-fragment Casida completed in {time.time() - self.time_begin:.2f}s")

        if use_coupling:
            self.coupled_results = self._compute_coupling()
            yield self.coupled_results
        else:
            yield self.gsystem.casida_results

    def _compute_coupling(self):
        """Inter-fragment coupling on rank 0; bcast coupled ``omega``, ``f``, ``Z`` to all ranks."""
        from casidapy.subsystem_coupling import run_subsystem_casida

        _, comm, is_mpi = _mpi_context(self.gsystem)
        nsub = len(self.drivers)

        xc_func = self.gsystem.total_evaluator.funcdicts.get("XC", None)
        ke_func = self.gsystem.total_evaluator.funcdicts.get("KE", None)
        if xc_func is None:
            sprint("ERROR: XC functional required for subsystem Casida coupling.")
            return None

        t0 = time.time()
        sprint("Computing subsystem Casida coupling...")

        sub_densities = _gather_subsystem_list(
            self.drivers,
            comm,
            is_mpi,
            nsub,
            lambda d: np.asarray(d.density),
        )
        serial_grid, serial_rho = _prepare_serial_coupling_fields(
            self.gsystem, comm, is_mpi,
        )

        coupled = None
        stream_paths = getattr(self, "_fragment_stream_paths", None)
        scratch_dir = getattr(self, "_casida_scratch_dir", None)
        if (not is_mpi) or comm.rank == 0:
            sprint("Coupling on rank 0 (serial global grid)...", flush=True)
            with _serial_coupling_gsystem(
                self.gsystem, serial_grid, serial_rho,
            ):
                coupled = run_subsystem_casida(
                    gsystem=self.gsystem,
                    drivers=self.drivers,
                    fragment_results=self.gsystem.casida_results,
                    xc_functional=xc_func,
                    ke_functional=ke_func,
                    sub_densities=sub_densities,
                    tda=self.options.get("tda", False),
                    rho_cutoff=self.options.get("rho_cutoff", 1e-3),
                    fxc_max=self.options.get("fxc_max", 20.0),
                    verbose=True,
                    fragment_stream_paths=stream_paths,
                )
        if is_mpi:
            comm.Barrier()

        disk_paths = getattr(self, "_fragment_disk_paths", None) or stream_paths
        cleanup_fragment_files(
            disk_paths if disk_paths else [],
            scratch_dir=scratch_dir,
            comm=comm,
            is_mpi=is_mpi,
        )
        if (not is_mpi) or comm.rank == 0:
            sprint("Casida fragment stream files removed.", flush=True)

        if (not is_mpi) or comm.rank == 0:
            excluded = getattr(self, "_uncoupled_excluded", None)
            if coupled is not None:
                coupled.update(
                    merge_coupled_and_uncoupled_spectrum(
                        coupled,
                        excluded,
                        fragment_results_full=getattr(
                            self, "_fragment_results_full", None,
                        ),
                        normalize_fosc=True,
                        tda=self.options.get("tda", False),
                    ),
                )

        if is_mpi:
            coupled = comm.bcast(_slim_coupled_results(coupled), root=0)

        sprint(f"Subsystem coupling completed in {time.time() - t0:.2f}s")

        if coupled is not None and coupled.get("omega_all") is not None:
            self.gsystem.casida_merged_spectrum = {
                "omega": coupled["omega_all"],
                "f": coupled["f_all"],
                "n_coupled": coupled["n_coupled"],
                "n_uncoupled": coupled["n_uncoupled"],
                "is_coupled": coupled["is_coupled"],
                "fragment_index": coupled["fragment_index"],
                "state_index": coupled["state_index"],
            }
        else:
            self.gsystem.casida_merged_spectrum = None

        self.gsystem.coupled_casida_results = coupled

        if (not is_mpi) or comm.rank == 0:
            excluded = getattr(self, "_uncoupled_excluded", None)
            if coupled is not None and "omega" in coupled:
                sprint(
                    f"Coupled excitation energies (eV):\n"
                    f"{coupled['omega'] * _HARTREE_TO_EV}",
                )
                f_cpl = np.asarray(coupled["f"], dtype=float)
                f_sum = float(np.nansum(f_cpl))
                sprint(
                    f"Coupled oscillator strengths (raw, sum={f_sum:.6g}):\n{f_cpl}",
                )
                if f_sum > 0.0 and np.isfinite(f_sum):
                    sprint(
                        "Coupled oscillator strengths (relative, sum=1):\n"
                        f"{f_cpl / f_sum}",
                    )
            if excluded is not None:
                _log_uncoupled_excluded_states(excluded)
                unc_txt = self.options.get(
                    "casida_uncoupled_txt", "casida_uncoupled_spectrum.txt",
                )
                if not os.path.isabs(unc_txt):
                    unc_txt = os.path.join(
                        os.environ.get("SLURM_SUBMIT_DIR") or os.getcwd(), unc_txt,
                    )
                write_uncoupled_excluded_txt(unc_txt, excluded)
                sprint(f"Uncoupled spectrum written to {unc_txt}", flush=True)
            if coupled is not None and coupled.get("omega_all") is not None:
                _log_merged_spectrum(coupled)
                png_path = _plot_merged_spectrum(coupled, self.options)
                if png_path and coupled is not None:
                    coupled["plot_path"] = png_path

        return coupled
