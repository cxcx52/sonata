"""Verify the 15 run directories and generate paper-ready tables from raw observations."""
import argparse
import csv
import hashlib
import json
from math import isclose, isfinite, prod
from pathlib import Path
from statistics import median

from check_sweep_wire import inspect, FilteredBall, BallIndex


ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def digest(path):
    return hashlib.sha256(path.read_text(encoding="utf-8").encode("utf-8")).hexdigest()


def require(test, message):
    if not test:
        raise ValueError(message)


def write_rows(path, rows):
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("runs", type=Path)
    parser.add_argument("sis", type=Path, help="Sage estimator output for the four optimized profiles")
    parser.add_argument("reference_sis", type=Path, help="Sage estimator output for the original profile")
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    out = args.output
    out.mkdir(exist_ok=True)
    files = ROOT/"parameters"
    profiles = read(files/"sweep-profiles.json")["profiles"]
    original = read(files/"ledger-two-point.json")
    sis = read(args.sis)
    old = read(args.reference_sis)
    require(sis["profile_ledger_lf_sha256"] == digest(files/"sweep-profiles.json"), "stale sweep estimates")
    require(sis["estimator_commit"] == old["estimator_commit"], "estimator revisions differ")
    machine = read(args.runs/"machine.json")
    topology = machine["topology"]
    cpus = [r["cpu"] for r in topology]
    require(len(cpus) == len(set(cpus)) == len({r["l3"] for r in topology}) == 3,
            "expected three trials on distinct CPUs and L3 caches")
    domain = FilteredBall(BallIndex(64,128),26)
    observations, groups, wires = [], [], []
    unique = set()
    for k,reference in [(16,False),(64,False),(256,False),(1024,False),(1024,True)]:
        name = "reference1024" if reference else f"k{k}"
        ledger = files/("ledger-two-point.json" if reference else "sweep-profiles.json")
        audits = old["results"] if reference else sis["profiles"][str(k)]["results"]
        ledger_roles = (profiles[str(k)]["representatives"] if not reference else
                        {name: [a,b,beta] for name,a,b,beta in original["roles"]
                         if name in original["embedding"].values()})
        require(len(ledger_roles) == 6 and ledger_roles.keys() == audits.keys(),
                f"{name}: incomplete representative set")
        for key,(rank,width,radius) in ledger_roles.items():
            r = audits[key]
            require((rank,width,radius) == (r["ring_rows"],r["ring_columns"],r["l2_radius"]),
                    f"{name}: stale {key} estimate")
        role, weak = min(audits.items(), key=lambda pair:pair[1]["log2_cost"]["core_classical"])
        hbits = (original["authenticated_node_union_bits"]["16384"] if reference
                 else profiles[str(k)]["issued_history_bits"])
        config = original if reference else profiles[str(k)]
        group = []
        for i,cpu in enumerate(cpus,1):
            directory = args.runs/f"{name}-r{i}"
            row = read(directory/"result.json")
            require(row["passed"] and row["k"] == k and row["reference"] == reference,
                    f"{directory.name}: wrong profile or verdict")
            require(row["affinity"] == [cpu] and row["binary_sha256"] == machine["binary_sha256"],
                    f"{directory.name}: affinity or executable changed")
            require(row["ledger_lf_sha256"] == digest(ledger), f"{directory.name}: wrong ledger")
            require(row["cpu"] == machine["cpu"] and row["kernel"] == machine["kernel"],
                    f"{directory.name}: hardware or kernel differs")
            measured = read(directory/"metrics.json")
            require(all(row[key] == value for key,value in measured.items()),
                    f"{directory.name}: result differs from C++ measurement")
            accounted = sum(row[key] for key in
                            ("verifier_coin_seconds","verifier_encode_seconds","verify_replay_seconds"))
            require(isclose(row["verifier_seconds"],accounted,rel_tol=1e-9) and
                    0 < row["prover_seconds"] < row["wall_seconds"] and
                    0 < row["verify_replay_seconds"] <= row["verifier_seconds"] < row["wall_seconds"]
                    and isfinite(row["wall_seconds"]) and row["max_rss_kib"] > 0,
                    f"{directory.name}: inconsistent process measurements")
            wire = inspect(directory, domain)
            for side in ("prover","verifier"):
                for suffix in ("bytes","sha256"):
                    key = side+"_"+suffix
                    require(wire[key] == row[key], f"{directory.name}: decoded {key} differs")
            require(wire["verifier_sha256"] not in unique,"reused verifier transcript")
            unique.add(wire["verifier_sha256"])
            wires.append(wire)
            obs = dict(profile=name, k=k, trial=i, cpu=cpu)
            for key in ("prover_seconds","verifier_seconds","verify_replay_seconds",
                        "verifier_coin_seconds","verifier_encode_seconds","setup_seconds",
                        "tamper_seconds","wall_seconds","prover_bytes","verifier_bytes","max_rss_kib"):
                obs[key] = row[key]
            obs.update(total_bytes=row["prover_bytes"]+row["verifier_bytes"],
                       seed_retries=sum(row["projection_attempts"])-5,
                       short_retries=sum(row["short_attempts"])-5)
            group.append(obs)
            observations.append(obs)
        summary = dict(profile=name, k=k, trials=len(group), weakest_role=role,
                       shared_crs_bytes=384*config["flat_CRS_cells"],
                       backend_pivot_call_coefficient=prod(2*(1+v["s"]) for v in config["layers"]),
                       classical_log2=weak["log2_cost"]["core_classical"],
                       quantum_sensitivity_log2=weak["log2_cost"]["core_quantum"],
                       history_statistical_bits=hbits)
        for key in ("prover_seconds","verifier_seconds","verify_replay_seconds","wall_seconds",
                    "prover_bytes","verifier_bytes","total_bytes","max_rss_kib"):
            values = [r[key] for r in group]
            for fn in (min,median,max):
                summary[key+"_"+fn.__name__] = fn(values)
        groups.append(summary)
    write_rows(out/"trials.csv",observations)
    write_rows(out/"summary.csv",groups)
    (out/"wire-checks.json").write_text(json.dumps(wires,indent=2)+"\n",encoding="utf-8")
    lines = ["# Five-layer ordinary transitions",
             "",
             "| Profile | Prover (s) | Verifier (s) | Public replay (s) | P→V (B) | V→P (B) | Total (B) | Peak RSS (MiB) | Classical / quantum model (log₂ work) |",
             "|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for r in groups:
        lines.append(f"| {r['profile']} | {r['prover_seconds_median']:.3f} | {r['verifier_seconds_median']:.3f} | "
                     f"{r['verify_replay_seconds_median']:.3f} | {r['prover_bytes_median']:,} | "
                     f"{r['verifier_bytes_median']:,} | {r['total_bytes_median']:,} | "
                     f"{r['max_rss_kib_max']/1024:.2f} | {r['classical_log2']:.3f} / "
                     f"{r['quantum_sensitivity_log2']:.3f} |")
    lines += ["", "Times and wire lengths are medians of three fresh interactions. "
              "Peak RSS is the largest individual-process high-water mark across those trials. "
              "Common public inputs and CRS are available to both parties before timing the prover. "
              "Verifier time includes issuing and encoding coins plus witness-free verification of the decoded files. "
              "Tamper checking is timed separately. Wire lengths count the two actual padded streams, including retries. "
              "Each profile ran as a wave of three trials pinned to separate physical cores and L3 caches on a shared host.",
              "", "The optimized profiles target 128-bit classical attack cost, with common screening floors of "
              "144 classical Core-SVP and 136 quantum Core-SVP sensitivity. The original 1024 control is reported "
              "under its original parameters. The theorem remains classical and conditional on the six timed Module-SIS "
              "advantages and the global extraction cutoff. The issued-attempt cap is 16,384."]
    (out/"experiment-table.md").write_text("\n".join(lines)+"\n", encoding="utf-8")
    print(f"Verified {len(unique)} distinct verifier transcripts across {len(groups)} profiles.")


if __name__ == "__main__":
    main()
