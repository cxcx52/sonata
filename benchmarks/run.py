"""Measure four reference kernels in separate Linux processes. Print JSON."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import random
from statistics import median
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]


def sample_challenges():
    sys.path.insert(0, str(ROOT))
    from extraction.ball_challenges import BallIndex, FilteredBall

    support = FilteredBall(BallIndex(64, 128), 26)
    rng = random.Random(5081)
    start = time.perf_counter_ns()
    rings = [support.sample(rng.getrandbits) for _ in range(21)]
    elapsed = time.perf_counter_ns() - start
    print(json.dumps(dict(wall_ns=elapsed, ring_count=len(rings), degree=64,
                          energy_cap=128, operator_cap=26,
                          sum_of_energies=sum(x*x for r in rings for x in r))))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path, help="built sonata_bench executable")
    parser.add_argument("--repeats", type=int, default=5)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("expected at least one sample per module")
    executable = args.executable.resolve(strict=True)
    if sys.platform != "linux" or not Path("/usr/bin/time").is_file():
        parser.error("requires Linux with GNU time installed in /usr/bin")

    cpuinfo = Path("/proc/cpuinfo").read_text().splitlines()
    cpu = next((line.split(":", 1)[1].strip() for line in cpuinfo
                if line.startswith("model name")), platform.machine())
    environment = dict(cpu=cpu, os=platform.freedesktop_os_release()["PRETTY_NAME"],
                       kernel=platform.release(), python=platform.python_version(),
                       affinity=sorted(os.sched_getaffinity(0)),
                       executable_sha256=hashlib.sha256(executable.read_bytes()).hexdigest())
    samples = {key: [] for key in ("encoding", "sumcheck", "projection", "challenge")}
    for _ in range(args.repeats):
        for name, rows in samples.items():
            cmd = ([sys.executable, __file__, "--challenge-worker"] if name == "challenge"
                   else [str(executable), name])
            run = subprocess.run(["/usr/bin/time", "-f", "max_rss_kib=%M", *cmd],
                                 capture_output=True, text=True)
            if run.returncode != 0:
                raise RuntimeError(f"{name} worker failed: {run.stderr.strip()}")
            row = json.loads(run.stdout)
            label, peak = run.stderr.strip().rsplit("=", 1)
            if label != "max_rss_kib":
                raise RuntimeError(f"unexpected GNU time result: {run.stderr!r}")
            rows.append(dict(row, peak_rss_kib=int(peak)))
    summary = {name: dict(median_ms=median(x["wall_ns"] for x in rows) / 1e6,
                          max_rss_kib=max(x["peak_rss_kib"] for x in rows))
               for name, rows in samples.items()}
    print(json.dumps(dict(environment=environment, summary=summary, samples=samples),
                     indent=2))


if __name__ == "__main__":
    if sys.argv[1:] == ["--challenge-worker"]:
        sample_challenges()
    else:
        main()
