import os as _os
import sys as _sys
import unittest as _unittest

from ._mpsym import *
from ._mpsym import __version__
from ._mpsym import __doc__

from . import _mpsym_tests


def _silent(cmd):
    _stdout = _sys.stdout
    _stderr = _sys.stderr

    try:
        with open(_os.devnull, 'w') as devnull:
            _sys.stdout = devnull
            _sys.stderr = devnull

            return cmd()

    finally:
        _sys.stdout = _stdout
        _sys.stderr = _stderr


def test(verbosity=-1):
    def run_tests(verbosity=0):
        suite = _unittest.TestLoader().loadTestsFromModule(_mpsym_tests)
        runner = _unittest.TextTestRunner(verbosity=verbosity)

        return runner.run(suite)

    if verbosity == -1:
        result = _silent(run_tests)
    else:
        result = run_tests(verbosity=verbosity)

    return 0 if result.wasSuccessful() else 1
