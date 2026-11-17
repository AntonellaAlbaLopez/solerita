"""solerita -- parametrize a ligand or a modified residue for AMBER, one
tutorial step at a time, seeing every command.

    from solerita import Parametrization
    p = Parametrization("CER", charge=0)
    p.inspect(); p.charges(); p.parameters(); p.library(); p.test()
    p.finish()

or, from a terminal:  solerita CER --charge 0
"""
__version__ = "1.0.0"

from .guide import Parametrization, StepFailed, amberhome

__all__ = ["Parametrization", "StepFailed", "amberhome", "__version__"]
