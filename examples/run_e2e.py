"""Run the five-layer example with reproducible short-challenge coins."""
import argparse
import os
from pathlib import Path
import random
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from extraction.ball_challenges import BallIndex, FilteredBall
from parameters.compute import OPERATOR_BOUND, PATH, RING_DEGREE


SUCCESS = "PASS: checkpoint verified through five layers and altered pivot rejected"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path, help="path to the built sonata_e2e")
    parser.add_argument("--run-full", action="store_true", required=True,
                        help="explicitly start the full-size computation")
    args = parser.parse_args()

    support = FilteredBall(BallIndex(degree=RING_DEGREE, energy=128), OPERATOR_BOUND)
    rng = random.Random(0x534F4E415441)
    requests = {(blocks, int(radix is None)) for blocks, radix in PATH}
    env = dict(os.environ, SONATA_E2E_CHALLENGE_PIPE="1")
    with subprocess.Popen([str(args.executable.resolve()), "--run-full"], env=env,
                          stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                          text=True, bufsize=1) as child:
        last_line = ""
        try:
            for line in child.stdout:
                last_line = line.rstrip("\n")
                if not line.startswith("SONATA_DRAW "):
                    print(line, end="", flush=True)
                    continue
                _, block_count, pivot_flag = line.split()
                blocks, pivot = int(block_count), int(pivot_flag)
                if (blocks, pivot) not in requests:
                    raise ValueError(f"invalid challenge request {line.strip()!r}")
                for j in range(blocks):
                    ring = support.sample(rng.getrandbits, nonzero=bool(pivot and j == 0))
                    child.stdin.write(" ".join(map(str, ring)) + "\n")
                child.stdin.flush()
            status = child.wait()
        finally:
            # Release the child even if sampling fails, the pipe breaks, or the
            # operator interrupts a large projection with Ctrl-C.
            if child.poll() is None:
                child.kill()
                child.wait()
    if status == 0 and last_line != SUCCESS:
        raise ValueError("executable exited without confirming checkpoint and pivot checks")
    return status if status >= 0 else 128 - status


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError) as exc:
        sys.exit(str(exc))
    except KeyboardInterrupt:
        sys.exit(130)
