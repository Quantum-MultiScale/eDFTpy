.. _ionization_potential:

===============================
Ionization Potential Tutorials
===============================

This tutorial computes the vertical ionization potential (IP) of a Cl⁻ ion embedded in a NaCl crystal,

.. math::

   \mathrm{IP} = E[\mathrm{Cl}] - E[\mathrm{Cl}^-],

by combining eDFTpy (embedding potential of the ion in the crystal) with PySCF (energies of the embedded ion with PBE, MP2 or CCSD).

The workflow has three steps:

1. **Standard run** (``standard/``): all ions in their ground state. Writes the embedding potential (EP) of one Cl⁻ and its screening potential.
2. **Oxidized run** (``oxidized/``): one Cl⁻ is replaced by a neutral Cl (the *impurity*). Writes the impurity-model EP of the oxidized ion.
3. **PySCF** (notebook ``ip_nacl.ipynb``): energies of Cl⁻ and Cl in vacuum, in the neutral EP and in the impurity EP, and the IP.

The embedding uses the Martyna–Tuckerman (MT) isolated boundary conditions for the subsystem's own Hartree and local pseudopotential, so that a charged subsystem does not interact with its own periodic images.

.. note::
   To avoid mixing files, run each step in its own folder, with the layout shown below.

Requirements
------------

* eDFTpy with Quantum ESPRESSO (QEpy) as the subsystem driver.
* PySCF with embedding-potential support (``vemb_mat``, ``mol.extemb``): https://github.com/JezsMartinez/pyscf (branch ``PRG_2026``).
* GBRV pseudopotentials ``na_pbe_v1.5.uspp.F.UPF`` and ``cl_pbe_v1.4.uspp.F.UPF`` in ``pseudo/``.

Folder layout
-------------

.. code-block:: text

   nacl_ip/
   ├── nacl.xyz
   ├── pseudo/            na_pbe_v1.5.uspp.F.UPF, cl_pbe_v1.4.uspp.F.UPF
   ├── common.py          PySCF driver (isolated + embedded passes)
   ├── ip_nacl.ipynb      step 3
   ├── standard/          input.ini, na.in, cl.in
   └── oxidized/          input.ini, na.in, cl.in, cl_oxidized.in

Input Files
-----------

The system is the rock-salt conventional cell (a = 5.64 Å, 4 Na⁺ + 4 Cl⁻, periodic).

.. dropdown:: nacl.xyz

   .. literalinclude:: /tutorials/advance_tutorial/ip_nacl/files/nacl.xyz

Quantum ESPRESSO templates for each subsystem (the charge of each ion is set by ``tot_charge`` and ``[MOL]``):

.. dropdown:: na.in

   .. literalinclude:: /tutorials/advance_tutorial/ip_nacl/files/na.in

.. dropdown:: cl.in

   .. literalinclude:: /tutorials/advance_tutorial/ip_nacl/files/cl.in

.. dropdown:: cl_oxidized.in

   .. literalinclude:: /tutorials/advance_tutorial/ip_nacl/files/cl_oxidized.in

Step 1: standard run
--------------------

.. dropdown:: standard/input.ini

   .. literalinclude:: /tutorials/advance_tutorial/ip_nacl/files/standard/input.ini

Each ion is its own subsystem (``decompose-method = distance``), giving 8 subsystems. The new keywords are:

* ``mt = True``: own Hartree and local pseudopotential with MT (isolated) boundary conditions.
* ``embedpot = .snpy``: write the embedding potential of the subsystem (``sub_cl_0.snpy``, ``sub_na_0.snpy``, …).
* ``mt_screening = True``: also write the screening potential ``sub_cl_0.screen.snpy``, used as the static reference in step 2.
* ``density-use_gaussians = True``: Gaussian core densities in the kinetic-energy functional.

Run it with one process per subsystem::

   $ cd standard
   $ mpirun -n 8 python -m edftpy input.ini

Step 2: oxidized run (impurity)
-------------------------------

.. dropdown:: oxidized/input.ini

   .. literalinclude:: /tutorials/advance_tutorial/ip_nacl/files/oxidized/input.ini

The first Cl⁻ (``cell-index = 4:5``) becomes its own subsystem ``[SUB_CL_c]``:

* ``density-ncharge = 0`` and ``density-nspin = 2``: neutral, open-shell Cl.
* ``mt_neutral_potential``: the screening potential of the same ion from step 1 (static reference of the impurity model).

Run it after step 1 has finished::

   $ cd oxidized
   $ mpirun -n 8 python -m edftpy input.ini

This writes ``sub_cl_c.snpy``, the impurity-model EP of the oxidized ion.

Step 3: IP with PySCF
---------------------

``common.py`` converts each ``.snpy`` EP to the PySCF grid, runs an isolated pass (no EP) and an embedded pass (EP as a one-electron potential), and returns:

* ``e_iso``: energy in vacuum,
* ``e_embd``: energy in the EP,
* ``e_fde``: interaction energy :math:`\mathrm{Tr}[D\,V_{\mathrm{emb}}]`,
* ``shift``: :math:`S = \frac{1}{N}\mathrm{Tr}[D\,(V^{\mathrm{static}} - V^{\mathrm{dynamic}})]`, with :math:`N` the number of electrons of the state.

.. dropdown:: common.py

   .. literalinclude:: /tutorials/advance_tutorial/ip_nacl/files/common.py
      :language: python

Open the notebook from ``nacl_ip/``, choose ``method`` (``pbe``, ``mp2`` or ``ccsd``) and run all cells:

.. toctree::
   :maxdepth: 1

   ip_nacl/ip_nacl

:download:`Download ip_nacl.ipynb <ip_nacl/ip_nacl.ipynb>` · :download:`Download common.py <ip_nacl/files/common.py>`

