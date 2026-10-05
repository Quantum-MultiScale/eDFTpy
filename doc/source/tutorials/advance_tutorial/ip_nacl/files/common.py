"""
Shared PySCF driver for the Ionic_Liquid IP/EA background jobs - ONE function
(run_state, used by run_and_save) that runs the isolated+embedded two-pass
scheme (same as the notebooks' pbe_iso/pbe) for PBE, MP2, or CCSD depending on
`method`. run_pbe.py / run_mp2.py / run_ccsd.py are thin CLI wrappers around
this, each launched as its own independent background process.
"""
import json
import os

import numpy as np
from pyscf import gto, dft, scf, cc
from dftpy.formats import io

#BASIS_BY_METHOD = {'pbe': 'def2-SVP', 'mp2': 'cc-pVDZ', 'ccsd': 'cc-pVDZ'}
#BASIS_BY_METHOD = {'pbe': 'def2-SVP', 'mp2': 'def2-SVP', 'ccsd': 'def2-SVP'}
BASIS_BY_METHOD = {'pbe': 'aug-cc-PVDZ', 'mp2': 'aug-cc-PVDZ', 'ccsd': 'aug-cc-PVDZ'}

GRID_LEVEL = 6


def _atomic_write_potential(final_path, data, ions):
    """io.write_potential(final_path, ...) via a per-process temp file +
    atomic rename - safe against multiple independent jobs regenerating the
    same shared multi-MB file at once (e.g. both IP-branch schemes read
    sub_ns2o4c2f6_0.snpy). The tmp marker goes BEFORE the extension, not
    after - dftpy's writer picks the file format from whatever follows the
    last dot, so a trailing ".tmp.<pid>" would be misread as the format
    itself."""
    root, ext = os.path.splitext(final_path)
    tmp_file = f"{root}.tmp{os.getpid()}{ext}"
    io.write_potential(tmp_file, data=data, ions=ions)
    os.replace(tmp_file, final_path)


def _to_pp(potential_file):
    """Auto-convert a .snpy embedding potential to .pp (what PySCF's spline
    reader wants), same as run_branch() in the IP_*.ipynb notebooks. Always
    regenerates from the .snpy rather than trusting a pre-existing .pp - an
    old .pp left over from an earlier run (different grid, different eDFTpy
    run) can silently disagree with the current .snpy otherwise, and
    dftpy's qepp reader has no way to detect that; it just fails (or worse,
    doesn't) trying to reshape stale data into the header's stale grid.

    Also writes a matching .xsf next to it (for visualizing the EP, e.g. in
    VESTA) as a side effect of this same conversion, so it's produced
    automatically by just running the notebook - no separate manual step."""
    if potential_file.endswith('.snpy'):
        pp_file = potential_file[:-len('.snpy')] + '.pp'
        xsf_file = potential_file[:-len('.snpy')] + '.xsf'
        ions, vemb_fft, _ = io.read(potential_file, kind='field')
        _atomic_write_potential(pp_file, vemb_fft, ions)
        _atomic_write_potential(xsf_file, vemb_fft, ions)
        return pp_file
    return potential_file


def _build_grid(geom_file, basis, spin, charge, grid_level=GRID_LEVEL):
    # spin/charge must be set BEFORE fromfile()/build() - fromfile() triggers an
    # internal validating build() immediately, using whatever spin/charge is
    # already on the mol at that point (default spin=0 would reject e.g. TFSI's
    # 137-electron neutral count or Pyr13's 73-electron one).
    mol = gto.M(basis=basis)
    mol.spin = spin
    mol.charge = charge
    mol.fromfile(geom_file)
    mol.build()
    grids = dft.gen_grid.Grids(mol)
    grids.level = grid_level
    grids.prune = None
    grids.build()
    return grids.coords, grids.weights


def _make_mf(mol, method):
    if method == 'pbe':
        mf = scf.UKS(mol)
        mf.xc = 'PBE'
    else:
        mf = scf.UHF(mol)
    return mf


def _run_isolated(method, geom_file, spin, charge, basis, extemb, coords, weights,
                   cell_cut=False, filename=None):
    """Isolated pass. For mp2/ccsd, vemb_mat() is read off the correlated
    object's own ._scf (matching uhf_mp2_ISO/uhf_ccsd_ISO in
    ~/Documents/1_IP/CCSDT_noneqe/{MP2_OS,CCSD_IS_py_f}.py) rather than off
    `mf` directly - vemb_mat() is a pure AO integral so this shouldn't change
    the matrix, but matching the validated reference exactly costs nothing."""
    mol = gto.M(basis=basis)
    mol.spin = spin
    mol.charge = charge
    mol.fromfile(geom_file)
    mol.extemb = extemb
    mol.ex_grids_coord = coords
    mol.ex_grids_weights = weights
    mol.build()

    mf = _make_mf(mol, method)
    mf.vemb = False
    mf.conv_tol = 1e-9
    mf.verbose = 0
    mf.max_cycle = 100
    mf.run()

    if method == 'mp2':
        corr = mf.MP2().run()
        e, scf_obj = corr.e_tot, corr._scf
    elif method == 'ccsd':
        corr = cc.UCCSD(mf).run()
        e, scf_obj = corr.e_tot, corr._scf
    else:
        e, scf_obj = mf.e_tot, mf

    mat, _ = scf_obj.vemb_mat(cell_cut=cell_cut, filename=filename)
    return e, mat


def _run_embedded(method, geom_file, spin, charge, basis, vemb_m):
    """Embedded pass. SCF convergence strategy and the density matrix used
    for the FDE interaction energy both depend on method, matching the
    validated reference scripts exactly (they are NOT the same settings as
    PBE's - no level_shift, DIIS on, looser conv_tol - and for mp2/ccsd the
    *correlated* relaxed density matrix is used, not the bare HF one)."""
    mol = gto.M(basis=basis)
    mol.spin = spin
    mol.charge = charge
    mol.fromfile(geom_file)
    mol.vemb_m = vemb_m
    mol.vemb = True
    mol.build()

    mf = _make_mf(mol, method)
    mf.conv_check = False
    mf.verbose = 2
    mf.max_cycle = 100
    mf.conv_tol = 1e-8
#    mf.level_shift = 0.0

    if method == 'pbe':
        mf = scf.addons.smearing_(mf, sigma=0.01, method="fermi")
        mf.run()
        e = mf.e_tot
        rdm1 = np.asarray(mf.make_rdm1())
    elif method == 'mp2':
        mf.run()
        corr = mf.MP2().run()
        e = corr.e_tot
        rdm1 = np.asarray(corr.make_rdm1(ao_repr=True))
    else:  # ccsd
        mf.run()
        corr = cc.UCCSD(mf).run()
        e = corr.e_tot
        rdm1 = np.asarray(corr.make_rdm1(ao_repr=True))

    fdee = np.einsum('imn,inm->', rdm1, vemb_m)
    return e, fdee, rdm1


def run_state(method, geom_file, spin, charge, extemb, grid_level=GRID_LEVEL,
              cell_cut=False, filename=None, shift_spec=None):
    """One charge/spin state, isolated + embedded passes, for the given method."""
    basis = BASIS_BY_METHOD[method]
    extemb = _to_pp(extemb)
    coords, weights = _build_grid(geom_file, basis, spin, charge, grid_level)
    e_iso, mat = _run_isolated(method, geom_file, spin, charge, basis, extemb,
                                coords, weights, cell_cut, filename)
    e_embd, fdee, rdm1 = _run_embedded(method, geom_file, spin, charge, basis, mat)
    out = dict(e_iso=e_iso, e_embd=e_embd, e_fde=fdee)
    if shift_spec is not None :
        out['shift'] = shift_correction(method, geom_file, spin, charge, rdm1, coords, weights, **shift_spec)
    return out


def run_and_save(method, out_file, geom_file,
                  standard_spin, standard_charge, perturbed_spin, perturbed_charge,
                  standard_extemb, perturbed_extemb, grid_level=GRID_LEVEL,
                  standard_cell_cut=False, standard_filename=None,
                  perturbed_cell_cut=False, perturbed_filename=None,
                  static_extemb=None, dynamic_extemb=None,
                  static_cell_cut=False, static_filename=None,
                  dynamic_cell_cut=False, dynamic_filename=None):
    """One IP or EA branch (standard + perturbed charge state), one method.
    Writes {out_file} on success - the notebook checks for this file to decide
    whether to skip re-running and just gather the result instead."""
    shift_spec = None
    if static_extemb is not None and dynamic_extemb is not None :
        shift_spec = dict(static_extemb=static_extemb, dynamic_extemb=dynamic_extemb,
                          static_cell_cut=static_cell_cut, static_filename=static_filename,
                          dynamic_cell_cut=dynamic_cell_cut, dynamic_filename=dynamic_filename)
    standard = run_state(method, geom_file, standard_spin, standard_charge,
                          standard_extemb, grid_level, standard_cell_cut, standard_filename, shift_spec)
    perturbed = run_state(method, geom_file, perturbed_spin, perturbed_charge,
                           perturbed_extemb, grid_level, perturbed_cell_cut, perturbed_filename, shift_spec)
    result = {'standard': standard, 'perturbed': perturbed}
    with open(out_file, 'w') as f:
        json.dump(result, f, indent=2)
    return result


def _vemb_matrix(method, geom_file, spin, charge, basis, extemb, coords, weights,
                 cell_cut=False, filename=None):
    """<mu|V|nu> of one embedding-potential file on the AO basis, shape (2, nao, nao)."""
    mol = gto.M(basis=basis)
    mol.spin = spin
    mol.charge = charge
    mol.fromfile(geom_file)
    mol.extemb = _to_pp(extemb)
    mol.ex_grids_coord = coords
    mol.ex_grids_weights = weights
    mol.build()
    mf = _make_mf(mol, method)
    mf.vemb = False
    mat, _ = mf.vemb_mat(cell_cut=cell_cut, filename=filename)
    return np.asarray(mat), mol.nelectron


def shift_correction(method, geom_file, spin, charge, rdm1, coords, weights,
                     static_extemb, dynamic_extemb,
                     static_cell_cut=False, static_filename=None,
                     dynamic_cell_cut=False, dynamic_filename=None):
    """S = 1/N * Tr[D (V_static - V_dynamic)]   (Ha)
    D: embedded rdm1 of this state; N = nelec of this state;
    V_static = neutral EP, V_dynamic = impurity EP."""
    basis = BASIS_BY_METHOD[method]
    v_static, nelec = _vemb_matrix(method, geom_file, spin, charge, basis, static_extemb,
                                   coords, weights, static_cell_cut, static_filename)
    v_dynamic, _ = _vemb_matrix(method, geom_file, spin, charge, basis, dynamic_extemb,
                                coords, weights, dynamic_cell_cut, dynamic_filename)
    return float(np.einsum('imn,inm->', rdm1, v_static - v_dynamic) / nelec)
