import numpy as np

from edftpy.mpi import sprint

def get_total_forces(drivers = None, gsystem = None, linearii=True, shift = True):
    # Global force: EWALD + local-PP + XC on the total density. This
    # value is REPLICATED on every MPI rank. Because the per-atom forces are
    # assembled with a final mp.vsum() over the global communicator, adding the
    # global force on every rank would count it nranks times. Add it on the
    # global root rank only, so it enters the reduction exactly once.
    global_forces = gsystem.get_forces(linearii = linearii)
    forces = np.zeros_like(global_forces)
    if gsystem.grid.mp.rank == 0 :
        forces[:] = global_forces
    # sprint('Total forces0 : \n', global_forces)
    for i, driver in enumerate(drivers):
        if driver is None : continue
        fs = driver.get_forces()
        ind = driver.subcell.ions_index
        if driver.technique == 'OF' :
            forces[ind] += fs
        elif driver.comm.rank == 0 :
            forces[ind] += fs
        # sprint('ind\n', ind, comm = driver.comm)
        # sprint('fs\n', fs, comm = driver.comm)
    forces = gsystem.grid.mp.vsum(forces)
    #-----------------------------------------------------------------------
    if shift :
        forces_shift = np.mean(forces, axis = 0)
        sprint('Forces shift :', forces_shift)
        forces -= forces_shift
    #-----------------------------------------------------------------------
    sprint('Total forces :')
    sprint(forces)
    return forces

def get_total_stress(drivers = None, gsystem = None, **kwargs):
    stress = gsystem.get_stress()
    for i, driver in enumerate(drivers):
        if driver is None : continue
        fs = driver.get_stress()
        # if driver.comm.rank > 0 : fs = 0.0
        # sprint('fs', fs, flush=True, comm=driver.comm)
        stress += fs
    stress = gsystem.grid.mp.vsum(stress)
    for i in range(2):
        for j in range(i+1, 3):
            stress[j,i] = stress[i,j]
    sprint('Total stress :')
    sprint(stress)
    return stress
