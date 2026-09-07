import numpy as np
import os

from dftpy.constants import ENERGY_CONV, STRESS_CONV
from edftpy.io import write, read
from edftpy.properties import get_electrostatic_potential

from edftpy.api.parse_config import config2optimizer, config2total_embed
from edftpy.mpi import graphtopo, sprint

from edftpy.utils.common import Field, Grid

import edftpy
import dftpy

def import_drivers(calcs = {}):
    """
    Import the engine of different drivers

    Notes:
        Must import driver firstly, before the mpi4py
    """
    if 'GSYSTEM' in calcs :
        config = calcs
        calcs = []
        for key in config :
            if key.startswith('SUB'):
                calc = config[key]['calculator']
                calcs.append(calc)
    fs = "{:>20s} Version : {}\n"
    info = fs.format('eDFTpy', edftpy.__version__)
    info += fs.format('DFTpy', dftpy.__version__)
    #
    try:
        from edftpy.engine import engine_qe
        if 'qe' in calcs or 'pwscf' in calcs or 'qepy' in calcs :
            info += fs.format('QEpy', engine_qe.__version__)
    except Exception as e:
        if 'qe' in calcs or 'pwscf' in calcs or 'qepy' in calcs : raise e
    try:
        from edftpy.engine import engine_castep
        if 'castep' in calcs :
            info += fs.format('Caspytep', engine_castep.__version__)
    except Exception as e:
        if 'castep' in calcs : raise e
    try:
        # Because Environ use same parallel technique as QE, need import after QEpy.
        from edftpy.engine import engine_environ
        if 'environ' in calcs :
            info += fs.format('Environ', engine_environ.__version__)
    except Exception as e:
        if 'environ' in calcs : return e
    return info

def conf2init(conf, parallel = False, **kwargs):
    info = import_drivers(conf)
    graphtopo = init_graphtopo(parallel, info = info, **kwargs)
    return graphtopo

def init_graphtopo(parallel = False, info = None, **kwargs):
    header = '*'*80 + '\n'
    if parallel :
        try :
            from mpi4py import MPI
            from mpi4py_fft import PFFT
            graphtopo.comm = MPI.COMM_WORLD
            header += 'Parallel version (MPI) on {0:>8d} processors\n'.format(graphtopo.comm.size)
        except Exception as e:
            raise e
    else :
        header += 'Serial version on {0:>8d} processor\n'.format(1)
    if graphtopo.is_root:
        #-----------------------------------------------------------------------
        # remove the stopfile
        if os.path.isfile('edftpy_stopfile'): os.remove('edftpy_stopfile')
        #-----------------------------------------------------------------------
    if info is None :
        fs = "{:>20s} Version : {}\n"
        header += fs.format('eDFTpy', edftpy.__version__)
        header += fs.format('DFTpy', dftpy.__version__)
    else :
        header += info
    header += '*'*80
    sprint(header)
    return graphtopo

def optimize_density_conf(config, **kwargs):
    opt = config2optimizer(config, **kwargs)
    opt.optimize()
    #-----------------------------------------------------------------------
    kefunc = opt.gsystem.total_evaluator.funcdicts.get('KE', None)
    if kefunc is not None and kefunc.name.startswith('MIX_'):
        opt.set_kedf_params(level = -1)
        for i, driver in enumerate(opt.drivers):
            if driver is None : continue
            driver.update_workspace(first = True)
        opt.optimize()
    #-----------------------------------------------------------------------
    energy = opt.energy
    sprint('Final energy (a.u.)', energy)
    sprint('Final energy (eV)', energy * ENERGY_CONV['Hartree']['eV'])
    sprint('Final energy (eV/atom)', energy * ENERGY_CONV['Hartree']['eV']/opt.gsystem.ions.nat)
    return opt

def optimize_embed(config, optimizer, lprint = False, mt = False, **kwargs):
    if not lprint :
        subkeys = [key for key in config if key.startswith('SUB')]
        for keysys in subkeys:
            if config[keysys].get("embedpot", None):
                lprint = True
                break
    if lprint :
        optimizer.set_kedf_params()
        # remove = ['PSEUDO', 'HARTREE', 'KE']
        remove = []
        removed = optimizer.gsystem.total_evaluator.update_functional(remove = remove)
        optimizer.set_global_potential()
        optimizer.gsystem.total_evaluator.update_functional(add = removed)
        global_nr = np.array(config["GSYSTEM"]["grid"]["nr"])
        has_cell_cut = False
        for subkey in config:
            if subkey.startswith('SUB') and "cell" in config[subkey] and "split" in config[subkey]["cell"]:
                cellsplit = config[subkey]["cell"]["split"]
                if cellsplit is not None and any(x > 0 for x in cellsplit):
                    has_cell_cut = True
                    break
        if not mt :
            subkeys = [key for key in config if key.startswith('SUB')]
            for keysys in subkeys:
                if  config[keysys].get("mt", None):
                    mt = True
                    break

        # mt_screening (generating the reference) and mt_neutral_potential (using it, impurity
        # model) both need a Hartree-only global_potential (global(periodic) restricted to the
        # Hartree functional alone) for every driver - see the per-driver loop below.
        # optimizer.gsystem.density is MPI-decomposed across the whole communicator, so this must
        # be computed once here via the existing, already-MPI-safe set_global_potential()
        # machinery - not inside the per-driver loop, where only some ranks execute and a
        # collective gather/broadcast would deadlock.
        need_hartree_only_global = mt and any(
            config[key].get("mt_screening", False) or config[key].get("mt_neutral_potential", None)
            for key in config if key.startswith('SUB')
        )
        hartree_only_global = {}
        if need_hartree_only_global :
            non_hartree_keys = [k for k in optimizer.gsystem.total_evaluator.funcdicts if k != 'HARTREE']
            removed_non_hartree = optimizer.gsystem.total_evaluator.update_functional(remove = non_hartree_keys)
            optimizer.set_global_potential()
            for i, other_driver in enumerate(optimizer.drivers):
                if other_driver is not None :
                    hartree_only_global[i] = other_driver.evaluator.global_potential.copy()
            optimizer.gsystem.total_evaluator.update_functional(add = removed_non_hartree)
            optimizer.set_global_potential()

        if has_cell_cut:
            global_embedding_potential = np.zeros(tuple(global_nr), dtype=float)
#            print("Global embedding potential initial", np.shape(global_embedding_potential))
            local_contributions = [(j, np.array(other_driver.evaluator.global_potential.data), tuple(other_driver.grid.nrR)) for j, other_driver in enumerate(optimizer.drivers) if other_driver is not None]
            all_contributions = graphtopo.comm.allgather(local_contributions)
            if graphtopo.is_root:
                for contributions in all_contributions:
                    for j, data, shape in contributions:
                        data_reshaped = data.reshape(shape)
                        index_j = optimizer.gsystem.graphtopo.graph.get_sub_index(j, in_global=True)
                        global_embedding_potential[index_j] = data_reshaped
            global_grid = Grid(optimizer.gsystem.grid.lattice, nr=tuple(global_nr))
            if graphtopo.is_root:
                global_embedding_potential = Field(global_grid, data=global_embedding_potential)
            else:
                global_embedding_potential = Field(global_grid)
            # Broadcast the data
            if graphtopo.is_root:
                data_to_broadcast = np.array(global_embedding_potential.data)
            else:
                data_to_broadcast = np.zeros(tuple(global_nr), dtype=float)
            graphtopo.comm.Bcast(data_to_broadcast, root=0)
            np.copyto(np.asarray(global_embedding_potential.data), data_to_broadcast)
#            print("Global embedding potential final", np.shape(global_embedding_potential))

        for i, driver in enumerate(optimizer.drivers):
            if driver is None : continue
            outfile = config[driver.key]["embedpot"]
            if outfile :
                driver = config2total_embed(config, driver = driver, optimizer = optimizer)
                if driver.technique == 'OF' or driver.comm.rank == 0 or graphtopo.isub is None:
                    removed_sub = driver.total_embed.update_functional(remove = remove)
                    index = optimizer.gsystem.graphtopo.graph.get_sub_index(i, in_global = True)
                    if has_cell_cut:
                        subsystem_potential = -driver.total_embed(driver.density, calcType = ['V']).potential
                        sprint("Using Cell-cut")
                        graph = optimizer.gsystem.graphtopo.graph

                        grid_shape = np.array(optimizer.gsystem.grid.nrR)
                        sub_shift = np.array(graph.sub_shift[i])
                        sub_shape = np.array(graph.sub_shape[i])

                        data_global = np.zeros(grid_shape)
                        data_local = np.array(subsystem_potential)
                        local_shape = data_local.shape
                        data_global[:local_shape[0],:local_shape[1],:local_shape[2]] = data_local
                        data = data_global
                        shape = data.shape

                        global_grid = Grid(optimizer.gsystem.grid.lattice, nr=tuple(grid_shape))
                        spacings = np.array(global_grid.spacings)

                        coords = np.indices(shape)
                        weights = np.abs(data)
                        current_center = np.array(local_shape) / 2.0
                        target_center = sub_shift + sub_shape / 2.0
                        shift_index = target_center - current_center

                        total_shift_real = shift_index * spacings

                        k = [2*np.pi*np.fft.fftfreq(shape[d], d=spacings[d]) for d in range(3)]
                        K = np.meshgrid(*k, indexing='ij')
                        phase = np.exp(-1j * (K[0]*total_shift_real[0] +
                                              K[1]*total_shift_real[1] +
                                              K[2]*total_shift_real[2] ))
                        data_fft = np.fft.fftn(data)
                        data = np.real(np.fft.ifftn(data_fft * phase))
                        subsystem_potential = Field(global_grid, data=data)

                        potential = global_embedding_potential - subsystem_potential
                        write(outfile, potential, optimizer.gsystem.ions, data_type = 'potential')
                    else:
                        if mt and config[driver.key].get("mt_neutral_potential", None):
                            fname = config[driver.key]["mt_neutral_potential"]
                            print("Reading MT neutral potential from file:", fname)
                            if graphtopo.is_root:
                                neutral_data = read(fname, kind='data', data_type='potential')
                            else:
                                neutral_data = None
                            if neutral_data is None:
                                subsystem_potential = driver.total_embed(driver.density, calcType=['V']).potential
                                potential = driver.evaluator.global_potential - subsystem_potential[index]
                                write(outfile, potential, driver.subcell.ions, data_type='potential')
                            else:
                                print("Using Impurity Model: ", fname)
                                grid_shape = np.array(optimizer.gsystem.grid.nrR)
                                global_grid = Grid(optimizer.gsystem.grid.lattice, nr=tuple(grid_shape))
                                # v_emb,imp = [v[n](KE+XC+PSEUDO+HARTREE) - v_bar[n_I](KE+XC+PSEUDO+HARTREE, MT)]
                                #           + own_hartree_mt - global_hartree_live + v_screen[n_neutral]
                                # own_hartree_mt cancels the Hartree already inside driver.total_embed's own
                                # potential (it includes Hartree, unlike the live evaluator), and
                                # global_hartree_live cancels the Hartree already inside global_potential,
                                # leaving v_emb,imp = v[n](KE+XC+PSEUDO) - v_bar[n_I](KE+XC+PSEUDO, MT) + v_screen[n_neutral].
                                neutral_screen_potential = Field(global_grid, data=neutral_data)
                                subsystem_potential = driver.total_embed(driver.density, calcType=['V']).potential
                                own_hartree_mt = driver.total_embed.funcdicts['HARTREE'](driver.density, calcType=['V']).potential
                                global_hartree_live = hartree_only_global[i]
                                potential = (driver.evaluator.global_potential - subsystem_potential[index]
                                             + own_hartree_mt[index]
                                             - global_hartree_live
                                             + neutral_screen_potential)
                                write(outfile, potential, driver.subcell.ions, data_type='potential')
                        else:
                            subsystem_potential = driver.total_embed(driver.density, calcType = ['V']).potential
                            potential = driver.evaluator.global_potential - subsystem_potential[index]
                            write(outfile, potential, driver.subcell.ions, data_type = 'potential')

                        if mt and config[driver.key].get("mt_screening", False) :
                            root, ext = os.path.splitext(outfile)
                            screen_outfile = root + ".screen" + ext
                            # v_screen[n] = v[n](periodic, global - all fragments) - v_bar[n_I](isolated/MT,
                            # own fragment) (Tolle et al. eq 11-12 via appendix eq A3): the Hartree-only
                            # slice of the ordinary embedding formula above, restricted to the Hartree
                            # functional alone. This is what mt_neutral_potential reads back in as the
                            # static neutral reference for the impurity model (EmbedEvaluator.get_embed_potential).
                            own_hartree = driver.total_embed.funcdicts['HARTREE'](driver.density, calcType=['V']).potential
                            screen_potential = hartree_only_global[i] - own_hartree[index]
                            write(screen_outfile, screen_potential, driver.subcell.ions, data_type='potential')
                    driver.total_embed.update_functional(add = removed_sub)
    return

def conf2output(config, optimizer):
    optimize_embed(config, optimizer)
    if config["GSYSTEM"]["density"]['output']:
        sprint("Write Density...")
        outfile = config["GSYSTEM"]["density"]['output']
        write(outfile, optimizer.density, optimizer.gsystem.ions)

    if config["OUTPUT"]["electrostatic_potential"]:
        sprint("Write electrostatic potential...")
        outfile = config["OUTPUT"]["electrostatic_potential"]
        v = get_electrostatic_potential(optimizer.gsystem)
        write(outfile, v, optimizer.gsystem.ions)

    for i, driver in enumerate(optimizer.drivers):
        if driver is None : continue
        outfile = config[driver.key]["density"]['output']
        if outfile :
            if driver.technique == 'OF' or driver.comm.rank == 0 or graphtopo.isub is None:
                write(outfile, driver.density, driver.subcell.ions)

    if "Force" in config["JOB"]["calctype"]:
        sprint("Calculate Force...")
        forces = optimizer.get_forces()
        ############################## Output Force ##############################
        if optimizer.gsystem.graphtopo.is_root:
            sprint("-" * 80)
            fabs = np.abs(forces)
            fmax, fmin, fave = fabs.max(axis = 0), fabs.min(axis = 0), fabs.mean(axis = 0)
            fmsd = (fabs * fabs).mean(axis = 0)
            fstr_f = " " * 4 + "{0:>12s} : {1:< 22.10f} {2:< 22.10f} {3:< 22.10f}"
            sprint(fstr_f.format("Max force (a.u.)", *fmax))
            sprint(fstr_f.format("Min force (a.u.)", *fmin))
            sprint(fstr_f.format("Ave force (a.u.)", *fave))
            sprint(fstr_f.format("MSD force (a.u.)", *fmsd))
            sprint("-" * 80)
        optimizer.gsystem.grid.mp.comm.Barrier()

    if "Stress" in config["JOB"]["calctype"]:
        sprint("Calculate Stress...")
        stress = optimizer.get_stress()
        ############################## Output Stress ##############################
        if optimizer.gsystem.graphtopo.is_root:
            sprint("-" * 80)
            sprint("Total stress (GPa):")
            stress_gpa = stress * STRESS_CONV["Ha/Bohr3"]["GPa"]
            sprint(stress_gpa)
            sprint("-" * 80)
        optimizer.gsystem.grid.mp.comm.Barrier()

    if config["OUTPUT"]["sub_temp"] :
        save = ['D', 'W']
    else :
        save = ['D']
    optimizer.stop_run(save = save)
    return

def __optimize_density_conf_test(config, **kwargs):
    opt = config2optimizer(config, **kwargs)
    opt.optimize()
    #-----------------------------------------------------------------------
    kefunc = opt.gsystem.total_evaluator.funcdicts.get('KE', None)
    if kefunc is not None and kefunc.name.startswith('MIX_'):
        opt.set_kedf_params(level = -1)
        for i, driver in enumerate(opt.drivers):
            if driver is None : continue
            driver.stop_run()
        opt2 = opt
        opt = config2optimizer(config, opt.gsystem.ions, graphtopo = opt.gsystem.graphtopo)
        #-----------------------------------------------------------------------
        for i, driver in enumerate(opt.drivers):
            if driver is None : continue
            if 'KE' in driver.evaluator.funcdicts :
                driver.evaluator.funcdicts['KE'].rhomax = opt2.drivers[i].evaluator.funcdicts['KE'].rhomax
            opt.gsystem.total_evaluator.funcdicts['KE'].rhomax = opt2.gsystem.total_evaluator.funcdicts['KE'].rhomax
        opt.optimize()
    #-----------------------------------------------------------------------
    energy = opt.energy
    sprint('Final energy (a.u.)', energy)
    sprint('Final energy (eV)', energy * ENERGY_CONV['Hartree']['eV'])
    sprint('Final energy (eV/atom)', energy * ENERGY_CONV['Hartree']['eV']/opt.gsystem.ions.nat)
    return opt
