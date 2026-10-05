.. _molecular_dynamics:

============================
Molecular Dynamics Tutorials
============================

This directory contains molecular dynamics tutorials using subsystem embedding. 
Source code and files for these examples can be found in the `eDFTpy GitHub Repository (Dynamics) <https://github.com/Quantum-MultiScale/eDFTpy/tree/dev/examples/dynamics>`_.

To run these examples, you need the eDFTpy configuration file (``input.ini``), an ASE Python script (``ase_nvt.py``) that uses the eDFTpy calculator, a SLURM submission script (``jobfile``, which depends on your local installation), the structure file (``*.xyz``), and pseudopotentials (all kept in the ``DATA/`` folder).

Nickel hydrated complex
========================
* **GitHub:** https://github.com/Quantum-MultiScale/eDFTpy/tree/dev/examples/dynamics/ni_h2o
* **System:** Nickel + water + 2 Cl⁻ ions (Coordinates: ``DATA/ni.xyz``, 342 atoms)
* **Type:** Molecular Dynamics (NVT) via ``ase_nvt.py``
* **Partition:** Distance-based decomposition (Ni: 0–19, Cl: 19–21, H₂O: 21–end)
* **Required files:** ``input.ini``, ``ase_nvt.py``, ``jobfile``, ``ni.xyz``, and pseudopotentials in ``DATA/``

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
   * - SUB_PD (Ni)
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

    sbatch jobfile


Palladium hydrated complex
==========================
* **GitHub:** https://github.com/Quantum-MultiScale/eDFTpy/tree/dev/examples/dynamics/pd_h2o
* **System:** Palladium + water + 2 Cl⁻ ions (Coordinates: ``DATA/pd.xyz``, 342 atoms)
* **Type:** Molecular Dynamics (NVT) via ``ase_nvt.py``
* **Partition:** Implicit distance-based partitioning (Pd: 0–13, Cl: 13–15, H₂O: 15–end)
* **Required files:** ``input.ini``, ``ase_nvt.py``, ``jobfile``, ``pd.xyz``, and pseudopotentials in ``DATA/``

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
   * - SUB_PD (Pd)
     - —
     - 40.0
     - 400.0
     - 0 : 13
   * - SUB_CL (Cl)
     - —
     - 40.0
     - 400.0
     - 13 : 15
   * - SUB_H2O
     - —
     - 40.0
     - 400.0
     - 15 :

**Run Command:**

.. code-block:: bash

    sbatch jobfile