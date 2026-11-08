"""Measure one full-profile run under Linux, including the Python challenge supplier."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import subprocess
import sys
import threading
import time


def rss_tree(pid):
    pending = [pid]
    seen = set()
    rss = 0
    while pending:
        current = pending.pop()
        if current in seen:
            continue
        seen.add(current)
        root = Path("/proc") / str(current)
        try:
            pending.extend(map(int, (root / f"task/{current}/children").read_text().split()))
            for line in (root / "status").read_text().splitlines():
                if line.startswith("VmRSS:"):
                    rss += int(line.split()[1])
                    break
        except (FileNotFoundError, ProcessLookupError):
            pass
    return rss


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("profile", choices=["window", "transition"])
    parser.add_argument("executable", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(exist_ok=False)
    executable = args.executable.resolve()
    cpu = next(line.split(":", 1)[1].strip() for line in Path("/proc/cpuinfo").read_text().splitlines()
               if line.startswith("model name"))
    cmd = [sys.executable, "-B", "examples/run_e2e.py", str(executable), "--run-full"]
    if args.profile == "transition":
        cmd.append("--transition")
    meta = dict(profile=args.profile, cpu=cpu, kernel=platform.release(), python=sys.version.split()[0],
                affinity=sorted(os.sched_getaffinity(0)), command=cmd,
                binary_sha256=hashlib.sha256(executable.read_bytes()).hexdigest(),
                start_utc=datetime.now(timezone.utc).isoformat())
    (args.output / "environment.json").write_text(json.dumps(meta, indent=2)+"\n")
    events = []
    peak = 0
    stop = threading.Event()
    start = time.monotonic()
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

    def monitor():
        nonlocal peak
        while not stop.is_set():
            peak = max(peak, rss_tree(proc.pid))
            stop.wait(0.1)

    worker = threading.Thread(target=monitor)
    worker.start()
    try:
        with (args.output / "trace.log").open("w", buffering=1) as log:
            for line in proc.stdout:
                seconds = time.monotonic()-start
                events.append(dict(seconds=seconds, message=line.rstrip()))
                print(f"{seconds:10.3f} {line.rstrip()}", file=log)
        code = proc.wait()
        elapsed = time.monotonic()-start
    finally:
        stop.set()
        worker.join()
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    phases = [e["message"] for e in events]
    verified = all(phases.count(f"{stage} backend layer {i}/5") == 1
                   for stage in ("Proving", "Verifying") for i in range(1, 6))
    passed = code == 0 and verified and phases[-1].startswith("PASS")
    result = dict(exit_code=code, passed=passed, wall_seconds=elapsed,
                  largest_child_maxrss_kib=usage.ru_maxrss, sampled_tree_peak_rss_kib=peak,
                  rss_sample_period_s=0.1, child_user_s=usage.ru_utime, child_sys_s=usage.ru_stime,
                  events=events)
    (args.output / "result.json").write_text(json.dumps(result, indent=2)+"\n")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
