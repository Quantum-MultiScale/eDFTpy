import os
import unittest

import numpy as np

from edftpy.utils.common import Grid, Ions
from edftpy.evaluator import EmbedEvaluator
from edftpy.api.parse_config import config2embed_evaluator

data_path = os.environ.get('EDFTPY_DATA_PATH')
if not data_path : data_path = 'DATA/'
if not os.path.exists(data_path) : data_path = '../DATA/'


class TestCellCutMTDisabled(unittest.TestCase):
    """config2embed_evaluator currently does not support mt=True together with
    cell-cut (a cut subsystem's own grid is a genuinely smaller box than
    GSYSTEM's - self-consistent MT there needs GSYSTEM's own grid, which this
    function does not have access to). Rather than build something incorrect,
    it warns and falls back to plain (non-MT) embedding for that subsystem.

    This is a regression test for that fallback, not a test of self-consistent
    cell-cut+MT support (which does not exist yet). If that combination is
    added back, this test should be replaced with one that exercises it
    directly instead of checking for the disabled fallback.
    """

    @classmethod
    def setUpClass(cls):
        cls.pp_o = os.path.join(data_path, 'O.pbe-rrkjus.UPF')
        cls.pp_h = os.path.join(data_path, 'H.pbe-rrkjus.UPF')
        if not (os.path.exists(cls.pp_o) and os.path.exists(cls.pp_h)):
            raise unittest.SkipTest('O/H.pbe-rrkjus.UPF not found under ' + data_path)

    def _base_config(self, cell_split):
        return {
            'GSYSTEM': {'kedf': {'kedf': 'TF'}, 'exc': {'xc': 'LDA'}},
            'MATH': {'linearie': True},
            'SUB_H2O': {
                'kedf': {'kedf': 'TF'},
                'embed': ['KE', 'XC', 'HARTREE', 'PSEUDO'],
                'exttype': None,
                'opt': {},
                'calculator': 'qe',
                'mt': True,
                'cell': {'split': cell_split},
            },
        }

    def _grid_and_ions(self, nr):
        grid = Grid(lattice=np.eye(3) * 11.735, nr=nr, full=False, direct=True)
        numbers = [8, 1, 1]
        positions = [[2.0, 2.2, 2.0], [2.5, 2.6, 2.0], [1.7, 2.6, 2.3]]
        ions = Ions(numbers=numbers, positions=positions, cell=grid.lattice, charges=[0, 0, 0])
        return grid, ions

    def test_cell_cut_disables_mt_with_warning(self):
        config = self._base_config(cell_split=[0.8, 0.8, 0.8])
        grid, ions = self._grid_and_ions((32, 32, 32))
        pplist = {'O': self.pp_o, 'H': self.pp_h}

        ev, exttype = config2embed_evaluator(config, 'SUB_H2O', ions, grid, pplist=pplist)
        self.assertIsInstance(ev, EmbedEvaluator)

        # mt=True was requested but cell-split is set, so PSEUDO/HARTREE must have
        # been built WITHOUT an mt screening kernel (plain periodic electrostatics).
        # LocalPP/Pseudo stores it privately as _mt; Hartree exposes it as .mt.
        pseudo = ev.funcdicts['PSEUDO']
        self.assertIsNone(getattr(pseudo, '_mt', 'MISSING'),
                           "cell-cut fragment's PSEUDO should fall back to mt=None")
        hartree = ev.funcdicts['HARTREE']
        self.assertIsNone(getattr(hartree, 'mt', 'MISSING'),
                           "cell-cut fragment's HARTREE should fall back to mt=None")

    def test_no_cell_cut_keeps_mt_enabled(self):
        # control: without cell-split, mt=True must actually take effect - confirms
        # the disable path above is specific to cell-cut, not mt being ignored outright.
        config = self._base_config(cell_split=None)
        grid, ions = self._grid_and_ions((32, 32, 32))
        pplist = {'O': self.pp_o, 'H': self.pp_h}

        ev, exttype = config2embed_evaluator(config, 'SUB_H2O', ions, grid, pplist=pplist)
        pseudo = ev.funcdicts['PSEUDO']
        self.assertIsNotNone(getattr(pseudo, '_mt', None),
                              "non-cell-cut fragment's PSEUDO should keep its mt screening")
        hartree = ev.funcdicts['HARTREE']
        self.assertIsNotNone(getattr(hartree, 'mt', None),
                              "non-cell-cut fragment's HARTREE should keep its mt screening")

    def test_cell_cut_embedding_potential_is_usable(self):
        # the fallback must still produce a working (non-crashing, finite) embedding
        # potential - this is the actual regression the cluster hit: cell-cut+mt
        # should degrade gracefully to plain embedding, not error out.
        config = self._base_config(cell_split=[0.8, 0.8, 0.8])
        grid, ions = self._grid_and_ions((32, 32, 32))
        pplist = {'O': self.pp_o, 'H': self.pp_h}

        ev, exttype = config2embed_evaluator(config, 'SUB_H2O', ions, grid, pplist=pplist)
        from edftpy.utils.common import Field
        rho = Field(grid, data=np.random.default_rng(0).random((32, 32, 32)) * 0.01)
        ev.global_potential = Field(grid, data=np.zeros((32, 32, 32)))
        ev.get_embed_potential(rho, with_ke=True)
        pot = np.asarray(ev.embed_potential)
        self.assertEqual(pot.shape, (32, 32, 32))
        self.assertTrue(np.isfinite(pot).all())


if __name__ == '__main__':
    unittest.main()
