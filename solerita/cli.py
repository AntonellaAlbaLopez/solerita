"""Command line:  solerita SME --charge 0

Runs every tutorial step in order, pausing before each one so you can read
what it is about to do. Use --no-pause to run straight through.
"""
import argparse
import shutil
import sys

from . import __version__
from .guide import Parametrization, StepFailed, amberhome


def check_install():
    print(f"solerita {__version__}")
    ok = True
    for tool in ("antechamber", "parmchk2", "prepgen", "tleap", "sqm",
                 "sander"):
        where = shutil.which(tool)
        ok &= bool(where)
        print(f"  {tool:<12} {where or 'NOT FOUND'}")
    if ok:
        print(f"  AMBERHOME    {amberhome()}")
        print("Everything is in place.")
        return 0
    print("Some AmberTools programs are missing. Activate the environment "
          "first:\n    conda activate solerita")
    return 1


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="solerita",
        description="Parametrize a ligand or a modified residue for AMBER, "
                    "following the tutorial step by step.",
        epilog="examples:\n"
               "  solerita CER --charge 0          ligand from the PDB\n"
               "  solerita SME --charge 0          modified residue\n"
               "  solerita my_lig.pdb --charge -1 --resname LIG\n"
               "  solerita --check                 is everything installed?",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("source", nargs="?",
                    help="three-letter PDB code (e.g. CER) or a structure file")
    ap.add_argument("--charge", type=int,
                    help="net charge of the molecule (required)")
    ap.add_argument("--resname", help="residue name to use (default: the code)")
    ap.add_argument("--kind", choices=["ligand", "residue"],
                    help="normally detected from the PDB entry")
    ap.add_argument("--offline", action="store_true",
                    help="use the copy of CER/SME shipped with solerita")
    ap.add_argument("--no-pause", action="store_true",
                    help="do not wait for Enter between steps")
    ap.add_argument("--check", action="store_true",
                    help="check that AmberTools is installed, then stop")
    a = ap.parse_args(argv)
    if a.check:
        return check_install()
    if not a.source:
        ap.print_help()
        return 1
    if a.charge is None:
        ap.error("please give the net charge, e.g. --charge 0")
    try:
        amberhome()
        p = Parametrization(a.source, charge=a.charge, resname=a.resname,
                            kind=a.kind, offline=a.offline,
                            pause=not a.no_pause)
        p.run_all()
    except (StepFailed, RuntimeError, ValueError) as exc:
        print(f"\nStopped: {exc}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("\nStopped.", file=sys.stderr)
        return 130
    return 0


if __name__ == "__main__":
    sys.exit(main())
