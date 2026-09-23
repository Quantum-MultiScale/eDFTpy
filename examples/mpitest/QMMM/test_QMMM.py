#!/usr/bin/env python3
import numpy as np
from edftpy.interface import conf2init
from edftpy.api.parse_config import config2optimizer
from edftpy.config import read_conf

def test_QMMM():
    '''Original TEST in 20 \AA cell size for testing reducing it to 5 \AA'''
    fname = 'qmmm.in'
    config = read_conf(fname)
    graphtopo = conf2init(config, parallel = True)
    optimizer = config2optimizer(config, graphtopo = graphtopo)
    assert len(optimizer.drivers) == 2
    optimizer.optimize()
    print('Energy: ',optimizer.energy)
#    graphtopo.assert_check(np.isclose(optimizer.energy, -17.34156, atol = 1E-3))
#    graphtopo.assert_check(np.isclose(optimizer.energy, -17.338300521208282, atol = 1E-3)) # 20 \AA
    graphtopo.assert_check(np.isclose(optimizer.energy, -12.164558943203993, atol = 1E-3)) # 5 \AA


if __name__ == "__main__":
    test_QMMM()
