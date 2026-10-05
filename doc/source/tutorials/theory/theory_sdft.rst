.. _theory_sdft:

=============
Subsystem DFT
=============

sDFT is a divide-and-conquer approach where the total system (e.g., a model of bulk water) is split into interacting fragments, called subsystems. The electron density of the total system is partitioned into subsystem electron densities, :math:`\{\rho_I(\mathbf{r})\}`, and written as:

.. math::

    \rho(\mathbf{r}) = \sum\limits_{I=1}^{N_S} \rho_{I}(\mathbf{r}),

where :math:`N_S` is the total number of subsystems.

The total electronic energy in sDFT is defined as a functional of all the subsystem densities:

.. math::

    E[\{\rho_I\}] = \sum\limits_I T_{s}[\rho_I] + V_{\mathrm{nuc}}[\rho] + E_H[\rho] + E_{\mathrm{xc}}[\rho] + T_{s}^{\mathrm{nad}}[\{\rho_I\}],

where :math:`V_{\mathrm{nuc}}[\rho]`, :math:`E_H[\rho]`, and :math:`E_{\mathrm{xc}}[\rho]` are defined as in standard Kohn-Sham DFT.

The electronic energy in sDFT can also be expressed as the sum of additive and non-additive energy terms:

.. math::

    E[\{\rho_I\}] = \sum\limits_I E_{I}[\rho_I] + E^{\mathrm{nad}}[\{\rho_I\}].

The additive part is calculated using the ground state Kohn-Sham-DFT functional computed with the external potential of subsystem :math:`I`, :math:`v_{ext}^I(\mathbf{r})`. The total external potential is given by the sum of all subsystem external potentials, :math:`v_\mathrm{ext}(\mathbf{r})=\sum_I^{N_S}v_\mathrm{ext}^I(\mathbf{r})`.

**Non-additive energy terms:**

* Coulomb non-additive energy:
  
  .. math::
      E_{Coul}^{nadd}[\{\rho_I\}] = \sum_{I, J\neq I}^{Ns} \int \rho_{I}(r) v_{ext}^J(\mathbf{r}) dr + \left(E_{H}[\rho] - \sum\limits_{I=1}^{N_S} E_{H}[\rho_{I}]\right)

* Exchange-correlation non-additive energy:
  
  .. math::
      E_{xc}^{nadd}[\{\rho_I\}] = E_{xc}[\rho] - \sum\limits_{I=1}^{N_S} E_{xc}[\rho_{I}]

* Kinetic non-additive energy:
  
  .. math::
      T_{s}^{nadd}[\{\rho_I\}] = T_{s}[\rho] - \sum\limits_{I=1}^{N_S} T_{s}[\rho_{I}]

The additive part of the kinetic energy is evaluated exactly from the subsystem KS orbitals. The non-additive kinetic energy :math:`T_{s}^{nadd}[\{\rho_I\}]` is evaluated approximately with pure density functionals.

The variational problem in sDFT is given by the set of KS-like equations with constrained electron density:

.. math::

    \left( -\frac{\nabla^2}{2} + v_{\mathrm{eff}}^{I}[\rho_I](\mathbf{r}) + v_{\mathrm{emb}}^{I}[\rho_{I},\rho](\mathbf{r}) \right) \psi_{i}^I(\mathbf{r}) = \epsilon_{i}^I \psi_{i}^I(\mathbf{r})

where the subsystem-specific effective KS potential contains the intra-subsystem contributions:

.. math::

    v_{\mathrm{eff}}^{I}[\rho_I](\mathbf{r}) = v_{\mathrm{ext}}^I(\mathbf{r}) + v_{{H}}[\rho_I](\mathbf{r}) + v_{\mathrm{xc}}[\rho_I](\mathbf{r}).

The interaction of the electrons of subsystem :math:`I` with the environment is represented by an embedding potential:

.. math::

    v_{\mathrm{emb}}^{I}[\rho_I,\rho](\mathbf{r}) = \sum\limits_{J,J \neq I} \left[ v_{\mathrm{ext}}^{J}(\mathbf{r}) + v_{{H}}[\rho_J](\mathbf{r}) \right] + v_{\mathrm{xc}}^{\mathrm{nad}}[\{\rho_{I}\}](\mathbf{r}) + v_{T_s}^{\mathrm{nad}}[\{\rho_{I}\}](\mathbf{r}).

The non-additive exchange-correlation potential and the non-additive kinetic potential are defined as:

.. math::

    v_{\mathrm{xc}}^{\mathrm{nad}}[\{\rho_{I}\}(\mathbf{r}) = \frac{\delta E_{xc}[\{\rho_I\}]}{\delta \rho_I(\mathbf{r})} = \frac{\delta E_{\mathrm{xc}}[\rho]}{\delta \rho(\mathbf{r})} - \frac{\delta E_{\mathrm{xc}}[\rho_I]}{\delta \rho_I(\mathbf{r})}.

.. math::

    v_{T_s}^{\mathrm{nad}}[\{\rho_{I}\}](\mathbf{r}) = \frac{\delta T_s[\{\rho_I\}]}{\delta \rho_I(\mathbf{r})} = \frac{\delta T_s[\rho]}{\delta \rho (\mathbf{r})} - \frac{\delta T_s[\rho_I]}{\delta \rho_I(\mathbf{r})}.

In practice, the set of subsystem KS equations can either be solved iteratively using a **Freeze-and-Thaw** procedure, or simultaneously for all subsystems with the total density being updated at every SCF cycle.