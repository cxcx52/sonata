"""Exercise the small interactive verifier with fresh operating-system coins."""
import argparse
from pathlib import Path
import secrets
import sys

from run_e2e import serve_verifier


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path, help="path to the built sonata_interactive")
    args = parser.parse_args()
    success = "PASS: small live verifier accepted the honest history and rejected all mutations"
    return serve_verifier(args.executable, "--small-live", {(3, 0), (3, 1)},
                          secrets.randbits, success)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError) as error:
        sys.exit(str(error))
    except KeyboardInterrupt:
        sys.exit(130)
