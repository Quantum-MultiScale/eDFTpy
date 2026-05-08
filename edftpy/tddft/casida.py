import numpy as np
import time

from edftpy.mpi import sprint
from edftpy.optimizer import Optimization

class CasidaTDDFT(Optimization):
    """
    The Casida TDDFT class for static Casida calculation

    Notes:
        Runs Casida calculation without time steps.
    """

    def __init__(self, drivers=None, gsystem = None, options=None, optimizer = None, **kwargs):
        super().__init__(drivers=drivers, gsystem = gsystem, options=options)
        default_options = {
            "olevel": 2,
            "sdft": 'sdft',
            "number_of_states": 10,
            "number_of_bands": 10,
        }
        self.options = default_options
        if isinstance(options, dict):
            self.options.update(options)
        self.optimizer = optimizer

    def initialization(self):
        self.optimizer.optimize()
        for driver in self.drivers:
            if driver is not None :
                driver.save(save = ['W', 'D'])
                driver.task = 'casida'
                sprint("Am I dying in update_workspace?")
                driver.update_workspace(first = True, options=self.options)

    def optimize(self, **kwargs):
        self.run(**kwargs)

    def run(self, **kwargs):
        for item in self.irun(**kwargs):
            pass
        return item

    def irun(self, restart = None, **kwargs):
        self.initialization()
        #-----------------------------------------------------------------------
        self.time_begin = time.time()
        #-----------------------------------------------------------------------
        self.casida_results = []
        for driver in self.drivers:
            if driver is not None and hasattr(driver, 'casida_results') and driver.casida_results is not None:
                self.casida_results.append(driver.casida_results)
            else:
                self.casida_results.append(None)
        yield self.casida_results