"""Measure one ordinary transition in the four-arity sweep on Linux."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import secrets
import sys
from time import monotonic

from run_e2e import serve_verifier


ROOT = Path(__file__).resolve().parents[1]
SUCCESS = "PASS: verified five layers from the encoded transition and rejected the altered pivot"


def run_profile(k, reference):
    if reference:
        path = ROOT / "parameters/ledger-two-point.json"
        layers = json.loads(path.read_text(encoding="utf-8"))["layers"]
        return [(v["s"], v["seed"]["bytes"]) for v in layers], 126647, 44409, path
    path = ROOT / "parameters/sweep-profiles.json"
    row = json.loads(path.read_text(encoding="utf-8"))["profiles"][str(k)]
    return ([(v["s"], v["seed_bytes"]) for v in row["layers"]],
            row["no_retry_bytes"]["prover"], row["no_retry_bytes"]["verifier"], path)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def checked_metrics(out, k, reference, layers, p_bytes, v_bytes):
    data = json.loads((out / "metrics.json").read_text())
    require(data["k"] == k and data["reference"] is reference, "binary profile mismatch")
    for field in ("prover_seconds", "verifier_seconds", "verifier_coin_seconds",
                  "verifier_encode_seconds", "verify_replay_seconds", "setup_seconds", "tamper_seconds"):
        require(math.isfinite(data[field]) and data[field] >= 0, f"invalid {field}")
    require(data["prover_seconds"] > 0 and data["verify_replay_seconds"] > 0, "missing main timings")
    p_attempts, c_attempts = data["projection_attempts"], data["short_attempts"]
    require(len(p_attempts) == len(c_attempts) == len(layers) == 5, "missing layer counters")
    for (blocks, seed_bytes), seeds, short in zip(layers, p_attempts, c_attempts):
        require(all(type(v) is int and 1 <= v <= 160 for v in (seeds, short)), "invalid retry count")
        p_bytes += seeds + short - 2
        v_bytes += (seeds - 1)*seed_bytes + (short - 1)*blocks*20
    for who, filename, expected in (("prover", "p-to-v.bin", p_bytes), ("verifier", "v-to-p.bin", v_bytes)):
        wire = (out / filename).read_bytes()
        require(len(wire) == data[who+"_bytes"] == expected,
                f"{filename}: file, measured count and retry-adjusted encoder budget disagree")
        data[who+"_sha256"] = hashlib.sha256(wire).hexdigest()
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("k", type=int, choices=(16, 64, 256, 1024))
    parser.add_argument("output", type=Path, help="new per-run directory")
    parser.add_argument("--reference", action="store_true", help="original 1024-input five-layer profile")
    args = parser.parse_args()
    if sys.platform != "linux" or (args.reference and args.k != 1024):
        parser.error("requires Linux; --reference is valid only at k=1024")
    import resource
    binary = args.executable.resolve(strict=True)
    out = args.output.resolve()
    layers, p_bytes, v_bytes, ledger = run_profile(args.k, args.reference)
    allowed = {(blocks, int(i == len(layers)-1)) for i,(blocks, _) in enumerate(layers)}
    argv = [str(args.k), str(out), "reference" if args.reference else "optimized"]
    out.mkdir(parents=True, exist_ok=False)
    start_utc = datetime.now(timezone.utc).isoformat()
    start = monotonic()
    with (out / "trace.log").open("x", buffering=1) as log:
        status = serve_verifier(binary, argv, allowed, secrets.randbits, SUCCESS, trace=log)
    wall = monotonic() - start
    if status:
        return status
    largest_rss = max(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
                      resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    result = checked_metrics(out, args.k, args.reference, layers, p_bytes, v_bytes)
    cpu = next(line.split(":", 1)[1].strip() for line in Path("/proc/cpuinfo").read_text().splitlines()
               if line.startswith("model name"))
    result.update(passed=True, start_utc=start_utc, wall_seconds=wall, max_rss_kib=largest_rss,
                  cpu=cpu, kernel=platform.release(), affinity=sorted(os.sched_getaffinity(0)),
                  python=platform.python_version(), binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                  ledger_lf_sha256=hashlib.sha256(ledger.read_text().encode("utf-8")).hexdigest())
    (out / "result.json").write_text(json.dumps(result, indent=2)+"\n")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (KeyError, OSError, ValueError) as error:
        sys.exit(str(error))
    except KeyboardInterrupt:
        sys.exit(130)
