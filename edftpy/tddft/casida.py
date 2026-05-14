import numpy as np
import time

from edftpy.mpi import sprint
from edftpy.optimizer import Optimization

class CasidaTDDFT(Optimization):
    """
    Subsystem Casida LR-TDDFT with inter-fragment coupling.

    Runs per-fragment Casida calculations via CasidaPy, then computes
    inter-fragment coupling through the non-additive XC (+ kinetic) kernel
    and solves the coupled eigenvalue problem for the full excitation spectrum.
    """

    def __init__(self, drivers=None, gsystem = None, options=None, optimizer = None, **kwargs):
        super().__init__(drivers=drivers, gsystem = gsystem, options=options)
        default_options = {
            "olevel": 2,
            "sdft": 'sdft',
            "number_of_states": 10,
            "number_of_bands": None,
            "coupling": True,
            "tda": False,
            "rho_cutoff": 1e-3,
            "fxc_max": 20.0,
        }
        self.options = default_options
        if isinstance(options, dict):
            self.options.update(options)
        self.optimizer = optimizer
        self.coupled_results = None

    def initialization(self):
        self.optimizer.optimize()
        for idx, driver in enumerate(self.drivers):
            if driver is not None :
                driver.save(save = ['W', 'D'])
                driver.task = 'casida'
                driver.update_workspace(first = True, options=self.options)

    def optimize(self, **kwargs):
        self.run(**kwargs)

    def run(self, **kwargs):
        for item in self.irun(**kwargs):
            pass
        return item

    def irun(self, restart = None, **kwargs):
        self.initialization()
        self.time_begin = time.time()

        self.gsystem.casida_results = []
        for idx, driver in enumerate(self.drivers):
            if driver is not None :
                driver.save(save = ['W', 'D'])
                driver.task = 'casida'
                if hasattr(driver, 'casida_results') and driver.casida_results is not None:
                    self.gsystem.casida_results.append(driver.casida_results)
                else:
                    self.gsystem.casida_results.append(None)
                if driver.casida_results is not None:
                    for key, value in driver.casida_results.items():
                        if key == 'rho_transition':
                            n_trans = len(value) if value is not None else 0
                            sprint(f"Subsystem {idx} Casida: {key}: {n_trans} transition densities stored")
                        else:
                            sprint(f"Subsystem {idx} Casida: {key}: \n", value)

        t_frag = time.time() - self.time_begin
        sprint(f"Per-fragment Casida completed in {t_frag:.2f}s")

        if self.options.get('coupling', True):
            self.coupled_results = self._compute_coupling()
            yield self.coupled_results
        else:
            yield self.gsystem.casida_results

    def _compute_coupling(self):
        """Compute inter-fragment coupling and solve coupled Casida equation."""
        from casidapy.subsystem_coupling import run_subsystem_casida

        xc_func = self.gsystem.total_evaluator.funcdicts.get('XC', None)
        ke_func = self.gsystem.total_evaluator.funcdicts.get('KE', None)

        if xc_func is None:
            sprint("WARNING: No XC functional found in gsystem.total_evaluator, "
                   "skipping coupling.")
            return None

        t0 = time.time()
        sprint("Computing subsystem Casida coupling...")

        coupled = run_subsystem_casida(
            gsystem=self.gsystem,
            drivers=self.drivers,
            fragment_results=self.gsystem.casida_results,
            xc_functional=xc_func,
            ke_functional=ke_func,
            tda=self.options.get('tda', False),
            rho_cutoff=self.options.get('rho_cutoff', 1e-3),
            fxc_max=self.options.get('fxc_max', 20.0),
            verbose=True,
        )

        t_coupling = time.time() - t0
        sprint(f"Subsystem coupling completed in {t_coupling:.2f}s")

        if coupled is not None and 'omega' in coupled:
            omega_eV = coupled['omega'] * 27.211386
            sprint(f"Coupled excitation energies (eV):\n{omega_eV}")
            sprint(f"Coupled oscillator strengths:\n{coupled['f']}")

        self.gsystem.coupled_casida_results = coupled
        return coupled