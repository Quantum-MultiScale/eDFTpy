.. _opt_density:

==============================
Density Optimization Tutorials
==============================

This directory contains density optimization tutorials using subsystem embedding (sDFT and QM/MM). 
Source code and files for these examples can be found in the `eDFTpy GitHub Repository (Opt_density) <https://github.com/Quantum-MultiScale/eDFTpy/tree/dev/examples/opt_density>`_.

To run these examples, you need the eDFTpy configuration file (``input.ini``), Quantum ESPRESSO input templates (``qe_in*.in`` files) that provide extra QE keywords for each subsystem, a SLURM submission script (``jobfile``, which depends on your local installation), the structure file (``*.xyz``), and pseudopotentials (all kept in the ``DATA/`` folder).

sDFT
====

Copper hydrated complex
-----------------------
* **GitHub:** https://github.com/Quantum-MultiScale/eDFTpy/tree/dev/examples/opt_density/sDFT/cu_h2o_water_mpi
* **System:** Copper Hexahidrated + water + 2 Cl⁻ (Coordinates: ``cu_hexa.xyz``, 342 atoms)
* **Partition:** Distance-based decomposition (Cu: 0–19, Cl: 19–21, H₂O: 21–end)
* **Required files:** ``input_cs.ini``, ``qe_in*.in``, ``jobfile``, ``cu_hexa.xyz``, and pseudopotentials in ``DATA/``

.. list-table:: Parameters
   :widths: 20 20 20 20 20
   :header-rows: 1

   * - Name
     - grid-ecut (Ry)
     - ecutwfc (Ry)
     - ecutrho (Ry)
     - cell-index
   * - **Global**
     - 2700
     - —
     - —
     - 0 : 342
   * - SUB_CU (Cu)
     - —
     - 40.0
     - 400.0
     - 0 : 19
   * - SUB_CL (Cl)
     - —
     - 40.0
     - 400.0
     - 19 : 21
   * - SUB_H2O
     - —
     - 40.0
     - 400.0
     - 21 :

**Run Command:**

.. code-block:: bash

    cd opt_density/sDFT/cu_h2o_water_mpi
    source ~env/bin/activate
    mpirun -n 126 python -m edftpy input_cs.ini > log


Glucose embedded in water
-------------------------
* **GitHub:** https://github.com/Quantum-MultiScale/eDFTpy/tree/dev/examples/opt_density/sDFT/glucose_water_mpi
* **System:** Glucose + water (Coordinates: ``opt.xyz``)
* **Partition:** Cell-index partitioning (KS1: 0–24, KS2: 24–end)
* **Required files:** ``input.ini``, ``qe_in.in``, ``jobfile``, ``opt.xyz``, and pseudopotentials in ``DATA/``

.. list-table:: Parameters
   :widths: 20 20 20 20 20
   :header-rows: 1

   * - Name
     - grid-ecut (Ry)
     - ecutwfc (Ry)
     - ecutrho (Ry)
     - cell-index
   * - **Global**
     - 2700
     - —
     - —
     - 0 : 342
   * - SUB_KS1
     - —
     - 40
     - 400
     - 0 : 24
   * - SUB_KS2
     - —
     - 40
     - 400
     - 24 :

**Run Command:**

.. code-block:: bash

    cd opt_density/sDFT/glucose_water_mpi
    source ~env/bin/activate
    mpirun -n 202 python -m edftpy input.ini > log

Bulk water (MT variant) 
-----------------------
* **GitHub:** https://github.com/Quantum-MultiScale/eDFTpy/tree/dev/examples/opt_density/sDFT/h2o_bulk_mpi_MT
* **System:** Water bulk (192 molecules, ``00.xyz``)
* **Partition:** Single subsystem with distance decomposition
* **Required files:** ``input.ini``, ``qe_in.in``, ``jobfile``, ``00.xyz``, and pseudopotentials in ``DATA/``

.. list-table:: Parameters
   :widths: 20 20 20 20 20
   :header-rows: 1

   * - Name
     - grid-ecut (Ry)
     - ecutwfc (Ry)
     - ecutrho (Ry)
     - cell-index
   * - **Global**
     - (grid-nr 150x150x150)
     - —
     - —
     - 0 : 576
   * - SUB_KS
     - —
     - 40
     - 400
     - 0 : 192

**Run Command (MT Variant):**

.. code-block:: bash

    cd opt_density/sDFT/h2o_bulk_mpi_MT
    source ~env/bin/activate
    mpirun -n 64 python -m edftpy input.ini > log


Sodium chloride in water
------------------------
* **GitHub:** https://github.com/Quantum-MultiScale/eDFTpy/tree/dev/examples/opt_density/sDFT/h2o_nacl_mpi
* **System:** Sodium + chloride + water (Coordinates: ``salt_h2o_2.xyz``)
* **Partition:** Cell-index partitioning (Na: 0–1, Cl: 1–2, H₂O: 2–end)
* **Required files:** ``input.ini``, ``jobfile``, ``salt_h2o_2.xyz``, and pseudopotentials in ``DATA/``

.. list-table:: Parameters
   :widths: 20 20 20 20 20
   :header-rows: 1

   * - Name
     - grid-ecut (Ry)
     - ecutwfc (Ry)
     - ecutrho (Ry)
     - cell-index
   * - **Global**
     - 1200
     - —
     - —
     - 0 : N
   * - SUB_NA (Na)
     - 2400
     - —
     - —
     - 0 : 1
   * - SUB_CL (Cl)
     - 2400
     - —
     - —
     - 1 : 2
   * - SUB_H2O
     - 2400
     - —
     - —
     - 2 :

**Run Command:**

.. code-block:: bash

    cd opt_density/sDFT/h2o_nacl_mpi
    source ~env/bin/activate
    mpirun -n 4 python -m edftpy input.ini | tee log


Palladium chloride in water
---------------------------
* **GitHub:** https://github.com/Quantum-MultiScale/eDFTpy/tree/dev/examples/opt_density/sDFT/pd_cl_water_mpi
* **System:** Palladium + chloride + water (Coordinates: ``pd.xyz``, 342 atoms)
* **Required files:** ``input.ini``, ``qe_in*.in``, ``jobfile``, ``pd.xyz``, and pseudopotentials in ``DATA/``

.. list-table:: Parameters
   :widths: 20 20 20 20 20
   :header-rows: 1

   * - Name
     - grid-ecut (Ry)
     - ecutwfc (Ry)
     - ecutrho (Ry)
     - cell-index
   * - **Global**
     - 2700
     - —
     - —
     - 0 : 342
   * - SUB_PD
     - —
     - 40
     - 400
     - 0 : 13
   * - SUB_CL
     - —
     - 40
     - 400
     - 13 : 15
   * - SUB_H2O
     - —
     - 40
     - 400
     - 15 :

**Run Command:**

.. code-block:: bash

    cd opt_density/sDFT/pd_cl_water_mpi
    source ~env/bin/activate
    mpirun -n 160 python -m edftpy input.ini > log


Water dimer
-----------
* **GitHub:** https://github.com/Quantum-MultiScale/eDFTpy/tree/dev/examples/opt_density/sDFT/h2o_dimer_mpi
* **System:** Water dimer (Coordinates: ``h2o_2.xyz``, 6 atoms)
* **Partition:** Single subsystem (0–6)
* **Required files:** ``edftpy_2.ini`` (adaptive) or ``edftpy_1.ini`` (non-adaptive), ``qe_in.in``, xyz coordinates, and pseudopotentials in ``DATA/``

.. list-table:: Parameters
   :widths: 20 20 20 20 20
   :header-rows: 1

   * - Name
     - grid-ecut (Ry)
     - ecutwfc (Ry)
     - ecutrho (Ry)
     - cell-index
   * - **Global**
     - 1600
     - —
     - —
     - 0 : 6
   * - SUB_HO
     - —
     - 40
     - 400
     - 0 : 6

**Adaptive Run:**

.. code-block:: bash

    cd opt_density/sDFT/h2o_dimer_mpi/adaptive
    source ~env/bin/activate
    mpirun -n 4 python -m edftpy edftpy_2.ini | tee log.2

**No-Adaptive Run:**

.. code-block:: bash

    cd opt_density/sDFT/h2o_dimer_mpi/no-adaptive
    source ~env/bin/activate
    mpirun -n 4 python -m edftpy edftpy_1.ini | tee log.1

**QE Reference Run:**

.. code-block:: bash

    cd opt_density/sDFT/h2o_dimer_mpi/qe
    source ~env/bin/activate
    mpirun -n 4 python -m qepy --pw.x -i qe_ho.in | tee log.qe



Water dimer (Jupyter notebook) 
------------------------------
* **GitHub:** https://github.com/Quantum-MultiScale/eDFTpy/tree/dev/examples/opt_density/sDFT/h2o_dimer_JN
* **Type:** Density Optimization (sDFT, single subsystem) via Jupyter notebook
* **Partition:** Single subsystem
* **Required files:** ``test_opt_one_subsystem.ipynb``, ``water.xsf``, and pseudopotentials in ``DATA/``

**Run Command:**

.. code-block:: bash

    cd opt_density/sDFT/h2o_dimer_JN
    jupyter notebook test_opt_one_subsystem.ipynb


QM/MM
=====

Glucose embedded in water
-------------------------
* **GitHub:** https://github.com/Quantum-MultiScale/eDFTpy/tree/dev/examples/opt_density/qmmm/glucose_water
* **System:** Glucose + water (Coordinates: ``511.xyz``)
* **Type:** Density Optimization (QM/MM) via ``edftpy qmmm.ini``
* **Partition:** Cell-index partitioning (KS1: 0–24, MM: 24–end)
* **Required files:** ``qmmm.ini``, ``qe_in.in``, ``jobfile``, ``511.xyz``, pseudopotentials, and MM parameters in ``DATA/``

.. list-table:: Parameters
   :widths: 20 20 20 20 20
   :header-rows: 1

   * - Name
     - grid-ecut (Ry)
     - ecutwfc (Ry)
     - ecutrho (Ry)
     - cell-index
   * - **Global**
     - 2700
     - —
     - —
     - 0 : 342
   * - SUB_KS1 (QM)
     - —
     - 40
     - 400
     - 0 : 24
   * - SUB_MM
     - —
     - —
     - —
     - 24 :

**Run Command:**

.. code-block:: bash

    cd opt_density/qmmm/glucose_water
    source ~env/bin/activate
    mpirun -n 64 python -m edftpy qmmm.ini > log


Water dimer
-----------
* **GitHub:** https://github.com/Quantum-MultiScale/eDFTpy/tree/dev/examples/opt_density/qmmm/h2o_dimer
* **System:** Water dimer (Coordinates: ``test.xyz``)
* **Type:** Density Optimization (QM/MM) via ``test_QMMM.py``
* **Partition:** Cell-index partitioning (KS1: 0–3, MM: 3–end)
* **Required files:** ``test.xyz``, and pseudopotentials in ``DATA/``

.. list-table:: Parameters
   :widths: 20 20 20 20 20
   :header-rows: 1

   * - Name
     - grid-ecut (Ry)
     - ecutwfc (Ry)
     - ecutrho (Ry)
     - cell-index
   * - **Global**
     - 2700
     - —
     - —
     - 0 : 6
   * - SUB_KS1 (QM)
     - —
     - 70
     - 400
     - 0 : 3
   * - SUB_MM
     - —
     - —
     - —
     - 3 :

**Run Command:**

.. code-block:: bash

    cd opt_density/qmmm/h2o_dimer/
    source ~env/bin/activate
    python test_QMMM.py | tee log.0


