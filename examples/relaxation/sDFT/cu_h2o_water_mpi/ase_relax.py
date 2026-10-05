#!/usr/bin/env python3
import os
import numpy as np

from edftpy.api.api4ase import eDFTpyCalculator
from edftpy.config import read_conf
from edftpy.interface import conf2init
from edftpy.mpi import pmi, sprint

import ase.io
from ase.io.trajectory import Trajectory
from ase.optimize import LBFGS, BFGS

np.random.seed(8888)

# Read eDFTpy configuration
conf = read_conf("./input_cs.ini")

# Initialize eDFTpy graph topology
graphtopo = conf2init(conf, pmi.size > 0)

# Read structure file from input.ini
cell_file = conf["PATH"]["cell"] + os.sep + conf["GSYSTEM"]["cell"]["file"]
atoms = ase.io.read(cell_file,format='extxyz')

# Attach eDFTpy calculator
calc = eDFTpyCalculator(config=conf, graphtopo=graphtopo)
#atoms.set_calculator(calc)
atoms.calc = calc

# Run relaxation
trajfile = "opt.traj"
#opt = LBFGS(atoms, trajectory=trajfile, memory=10, use_line_search=False)
opt = BFGS(atoms,
        trajectory = trajfile, logfile = 'opt.log')

# Relax 
#opt.run(fmax=0.3)
opt.run(fmax = 0.02,steps=400)

# Save final relaxed structure
traj = Trajectory(trajfile)
ase.io.write("opt.xyz", traj[-1])

sprint("Relaxation finished. Final structure written to opt.xyz")
