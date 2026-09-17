.. _theory_qmmm:

====================================
A Density-Functionalization of QM/MM
====================================

The central idea is to assign an electron density to both the QM subsystem, :math:`\rho_{QM}(\mathbf{r})`, and the MM subsystem, :math:`\rho_{MM}(\mathbf{r})`, with the total electron density given by their sum and the energy functional borrowed from rigorous sDFT:

.. math::

    \rho(\mathbf{r}) = \rho_{QM}(\mathbf{r}) + \rho_{MM}(\mathbf{r})

.. math::

    E[\rho_{QM},\rho_{MM}] = E[\rho_{QM}] + E[\rho_{MM}] + E^{nad}[\rho_{QM},\rho_{MM}].

The additive part of the energy is given by:

.. math::

    E[\rho_{QM}] = T_s[\rho_{QM}] + E_H[\rho_{QM}] + E_{xc}[\rho_{QM}] + \int v_{QM}(\mathbf{r})\rho_{QM}(\mathbf{r}) d\mathbf{r}

.. math::

    E[\rho_{MM}] = T_s[\rho_{MM}] + E_H[\rho_{MM}] + E_{xc}[\rho_{MM}] + \int v_{MM}(\mathbf{r})\rho_{MM}(\mathbf{r}) d\mathbf{r}

The nonadditive energy is thus given by:

.. math::

    E^{nad}[\rho_{QM},\rho_{MM}] = T_s^{nad}[\rho_{QM},\rho_{MM}] + E_{xc}^{nad}[\rho_{QM},\rho_{MM}] + \int \frac{\rho_{QM}(\mathbf{r})\rho_{MM}(\mathbf{r})}{|\mathbf{r}-\mathbf{r}'|}d\mathbf{r} d\mathbf{r}' + \int \rho_{MM}(\mathbf{r})v_{QM}(\mathbf{r})d\mathbf{r} + \int \rho_{QM}(\mathbf{r})v_{MM}(\mathbf{r})d\mathbf{r}

For each MM site :math:`i`, we represent the valence electron density by a Gaussian function centered at the site position, :math:`\mathbf{R}_i`, with an adjustable width :math:`\sigma_i`:

.. math::

    \rho_{q_i}(\mathbf{r}) = (N_i - q_i) g_{\sigma_i}(\mathbf{r} - \mathbf{R}_i)

The total valence electron density is then given by the sum of these contributions:

.. math::

    \rho_{MM}(\mathbf{r}) = \sum_{i\in \text{MM charges}} \rho_{q_i}(\mathbf{r}) + \sum_{j\in \text{MM dipoles}} \rho_{\mu_j}(\mathbf{r}).

**QM embedding potential:**

.. math::

    v_{emb}^{QM}(\mathbf{r}) = v_{MM}(\mathbf{r}) + \int d\mathbf{r}' \frac{\rho_{MM}(\mathbf{r}')}{|\mathbf{r}-\mathbf{r}'|} + \frac{\delta T_s^{nad}}{\delta \rho_{QM}(\mathbf{r})} + \frac{\delta E_{xc}^{nad}}{\delta \rho_{QM}(\mathbf{r})}

**MM embedding potential:**

.. math::

    v_{emb}^{MM}(\mathbf{r}) = v_{QM}(\mathbf{r}) + \int d\mathbf{r}' \frac{\rho_{QM}(\mathbf{r}')}{|\mathbf{r}-\mathbf{r}'|} + \frac{\delta T_s^{nad}}{\delta \rho_{MM}(\mathbf{r})} + \frac{\delta E_{xc}^{nad}}{\delta \rho_{MM}(\mathbf{r})}