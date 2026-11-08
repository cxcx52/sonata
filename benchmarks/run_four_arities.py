"""Run the four optimized arities in identical three-core waves."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import signal
import subprocess
import sys
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("output", type=Path, help="new results directory")
    parser.add_argument("--cpus", nargs=3, required=True, type=int)
    args = parser.parse_args()
    if len(set(args.cpus)) != 3 or not set(args.cpus) <= os.sched_getaffinity(0):
        parser.error("choose three distinct available logical CPUs")
    topology = []
    for cpu in args.cpus:
        node = Path(f"/sys/devices/system/cpu/cpu{cpu}")
        def read(relative):
            return (node/relative).read_text().strip()
        topology.append(dict(cpu=cpu, socket=read("topology/physical_package_id"),
                             core=read("topology/core_id"), siblings=read("topology/thread_siblings_list"),
                             l3=read("cache/index3/shared_cpu_list")))
    if len({item["l3"] for item in topology}) != 3:
        parser.error("select three CPUs with distinct last-level caches")
    binary = args.binary.resolve(strict=True)
    out = args.output.resolve()
    out.mkdir()
    cpu_name = next(x.split(":",1)[1].strip() for x in Path("/proc/cpuinfo").read_text().splitlines()
                    if x.startswith("model name"))
    manifest = dict(start_utc=datetime.now(timezone.utc).isoformat(), cpu=cpu_name,
                    architecture=platform.machine(), kernel=platform.release(), python=platform.python_version(),
                    topology=topology, binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                    schedule="3 simultaneous trials per wave; wait for all before next wave")
    (out/"machine.json").write_text(json.dumps(manifest, indent=2)+"\n")
    for k in (16,64,256,1024):
        name = f"k{k}"
        children = []
        start = time.monotonic()
        print(f"START {name} {datetime.now(timezone.utc).isoformat()}",flush=True)
        try:
            for i,cpu in enumerate(args.cpus,1):
                path = out / f"{name}-r{i}"
                cmd = ["taskset","-c",str(cpu),sys.executable,"-B","examples/run_sweep.py",
                       str(binary),str(k),str(path)]
                with (out/f"{path.name}.log").open("w") as log:
                    children.append(subprocess.Popen(cmd,stdin=subprocess.DEVNULL,stdout=log,
                                                     stderr=subprocess.STDOUT, start_new_session=True))
            codes = [child.wait() for child in children]
            if any(codes):
                raise SystemExit(f"{name}: exits {codes}")
        finally:
            for child in children:
                if child.poll() is None:
                    try:
                        os.killpg(child.pid, signal.SIGTERM)
                    except ProcessLookupError:
                        pass
                    child.wait()
        print(f"PASS {name} {time.monotonic()-start:.1f}s",flush=True)
    print("COMPLETE: 12 optimized trials",flush=True)


if __name__ == "__main__":
    main()
