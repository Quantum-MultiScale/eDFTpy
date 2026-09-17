.. _relaxation:

====================
Relaxation Tutorials
====================

This directory contains structural relaxation tutorials using subsystem embedding. 
Source code and files for these examples can be found in the `eDFTpy GitHub Repository (Relaxation) <https://github.com/Quantum-MultiScale/eDFTpy/tree/dev/examples/relaxation>`_.

To run these examples, you need the eDFTpy configuration file (``input.ini``), an ASE Python script (``ase_relax.py``) that uses the eDFTpy calculator, a SLURM submission script (``jobfile``, which depends on your local installation), the structure file (``*.xyz``), and pseudopotentials (all kept in the ``DATA/`` folder).

sDFT
====

Copper hydrated complex
-----------------------
* **GitHub:** https://github.com/Quantum-MultiScale/eDFTpy/tree/dev/examples/relaxation/sDFT/cu_h2o_water_mpi
* **System:** Copper Hexahidrated + water + 2 Cl⁻ (Coordinates: ``cu_hexa.xyz``, 342 atoms)
* **Type:** Relaxation (sDFT) via ``ase_relax.py``
* **Partition:** Distance-based decomposition (Cu: 0–19, Cl: 19–21, H₂O: 21–end)
* **Required files:** ``input.ini``, ``ase_relax.py``, ``jobfile``, ``cu_hexa.xyz``, and pseudopotentials in ``DATA/``

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
     - 40
     - 400
     - 19 : 21
   * - SUB_H2O
     - —
     - 40
     - 400
     - 21 :

**Run Command:**

.. code-block:: bash

    cd relaxation/sDFT/cu_h2o_water_mpi
    source ~env/bin/activate
    mpirun -n $SLURM_NTASKS python -m ase_relax.py > log


Copper hydrated complex (no adaptive dynamic adjustment)
--------------------------------------------------------
* **GitHub:** https://github.com/Quantum-MultiScale/eDFTpy/tree/dev/examples/relaxation/sDFT/cu_h2o_water_mpi_no_ada
* **System:** Copper Hexahidrated + water + 2 Cl⁻ (Coordinates: ``cu_hexa.xyz``, 342 atoms)
* **Type:** Relaxation (sDFT, no adaptive density adjustment)
* **Partition:** Distance-based (without ADA)
* **Required files:** ``input.ini``, ``ase_relax.py``, ``jobfile``, ``cu_hexa.xyz``, and pseudopotentials in ``DATA/``

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
     - 40
     - 400
     - 19 : 21
   * - SUB_H2O
     - —
     - 40
     - 400
     - 21 :

**Run Command:**

.. code-block:: bash

    cd relaxation/sDFT/cu_h2o_water_mpi_no_ada
    source ~env/bin/activate
    mpirun -n $SLURM_NTASKS python -m ase_relax.py > log


Glucose embedded in water
-------------------------
* **GitHub:** https://github.com/Quantum-MultiScale/eDFTpy/tree/dev/examples/relaxation/sDFT/glucose_water
* **System:** Glucose + water (Coordinates: ``pert_1.xyz``)
* **Type:** Relaxation (sDFT)
* **Partition:** Cell-index partitioning (KS1: 0–24, KS2: 24–end)
* **Required files:** ``input.ini``, ``ase_relax.py``, ``jobfile``, ``pert_1.xyz``, and pseudopotentials in ``DATA/``

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
     - 0 : N
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

    cd relaxation/sDFT/glucose_water
    module use /projects/community-old/modulefiles
    module load gcc/15.2.0/openssl/3.6.0-krasting intel/19.1.1 gcc/14.2.0-cermak libffi/3.3-gc563
    source /projectsn/mp1009_1/for_all/software_2026/edftpy_2026_2/bin/activate
    mpirun -n $SLURM_NTASKS python -m ase_relax.py > log


Palladium hydrated complex
--------------------------
* **GitHub:** https://github.com/Quantum-MultiScale/eDFTpy/tree/dev/examples/relaxation/sDFT/pd_complex
* **System:** Palladium + water + Cl⁻ (Coordinates: ``pd.xyz``, 342 atoms)
* **Type:** Relaxation (sDFT)
* **Partition:** Distance-based decomposition (Pd: 0–13, Cl: 13–15, H₂O: 15–end)
* **Required files:** ``input.ini``, ``ase_relax.py``, ``jobfile``, ``pd.xyz``, and pseudopotentials in ``DATA/``

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
     - 40.0
     - 400.0
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

    cd relaxation/sDFT/pd_complex
    source ~env/bin/activate
    mpirun -n $SLURM_NTASKS python -m ase_relax.py > log


Water cluster (small)
---------------------
* **GitHub:** https://github.com/Quantum-MultiScale/eDFTpy/tree/dev/examples/relaxation/sDFT/water_small
* **System:** Small water cluster (27 molecules, ``water_subcell.xyz``)
* **Type:** Relaxation (sDFT)
* **Partition:** Cell-index partitioning (KS: 0–24, KS8: 24–27)
* **Required files:** ``input.ini``, ``ase_relax.py``, ``jobfile``, ``water_subcell.xyz``, and pseudopotentials in ``DATA/``

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
     - 0 : 81
   * - SUB_KS
     - —
     - 40
     - 400
     - 0 : 24
   * - SUB_KS8
     - —
     - 40
     - 400
     - 24 : 27

**Run Command:**

.. code-block:: bash

    cd relaxation/sDFT/water_small
    source ~env/bin/activate
    mpirun -n $SLURM_NTASKS python -m ase_relax.py > log



QM/MM
=====
QM/MM embedding with pwscf for the QM region and mbx for the MM region.

Glucose embedded in water
-------------------------
* **GitHub:** https://github.com/Quantum-MultiScale/eDFTpy/tree/dev/examples/relaxation/qmmm/glucose_water
* **System:** Glucose + water (Coordinates: ``pert_1.xyz``, Species: C, O, H)
* **Type:** Relaxation (QM/MM) via ``ase_relax.py``
* **Partition:** Cell-index partitioning (QM atoms 0–24, MM atoms 24–end)
* **Required files:** ``input.ini``, ``ase_relax.py``, ``jobfile``, ``pert_1.xyz``, and pseudopotentials in ``DATA/``

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
     - 0 : N
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

    cd relaxation/qmmm/glucose_water
    source ~env/bin/activate
    mpirun -n $SLURM_NTASKS python -m ase_relax.py > log