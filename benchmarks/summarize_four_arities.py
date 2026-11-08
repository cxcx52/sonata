"""Verify twelve five-layer trials and tabulate the optimized configurations."""
import argparse
import csv
import hashlib
import json
from fractions import Fraction
from math import isclose, isfinite, log2, prod
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
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    out = args.output
    out.mkdir(exist_ok=True)
    files = ROOT/"parameters"
    ledger = files/"sweep-profiles.json"
    profiles = read(ledger)["profiles"]
    sis = read(args.sis)
    require(sis["profile_ledger_lf_sha256"] == digest(files/"sweep-profiles.json"), "stale sweep estimates")
    require(sis["estimator_commit"] == "53da5982597709ba0fdf94ea37a84d822310fd84", "estimator revision differs")
    machine = read(args.runs/"machine.json")
    topology = machine["topology"]
    cpus = [r["cpu"] for r in topology]
    require(len(cpus) == len(set(cpus)) == len({r["l3"] for r in topology}) == 3,
            "expected three trials on distinct CPUs and L3 caches")
    domain = FilteredBall(BallIndex(64,128),26)
    observations, groups, wires, attacks = [], [], [], []
    unique = set()
    for k in (16,64,256,1024):
        name = f"k{k}"
        config = profiles[str(k)]
        audits = sis["profiles"][str(k)]["results"]
        ledger_roles = config["representatives"]
        require(len(ledger_roles) == 6 and ledger_roles.keys() == audits.keys(),
                f"{name}: incomplete representative set")
        for key,(rank,width,radius) in ledger_roles.items():
            r = audits[key]
            require((rank,width,radius) == (r["ring_rows"],r["ring_columns"],r["l2_radius"]),
                    f"{name}: stale {key} estimate")
            block = r["bkz_block_size"]
            require(r["matzov_within_fit"] == (block <= 1024), f"{name}: {key} MATZOV fit flag")
            crt = r["single_CRT_component"]
            crt_block = crt.get("bkz_block_size")
            attacks.append(dict(k=k, role=key, ring_rows=rank, ring_columns=width, l2_radius=radius,
                                bkz_block=block, **r["log2_cost"],
                                matzov_fit="within range" if block <= 1024 else "extrapolated (block > 1024)",
                                crt_screen=crt["status"], crt_block=crt_block,
                                crt_block_in_matzov_range=None if crt_block is None else crt_block <= 1024))
        role, weak = min(audits.items(), key=lambda pair:pair[1]["log2_cost"]["core_classical"])
        history = 16384 * Fraction(config["zero_pending_exact"])
        conditional_budget = history + Fraction(1,2**130) + Fraction(6,2**132)
        require(history <= Fraction(1,2**130), f"{name}: issued-attempt statistical budget exceeded")
        # The CRS stores column-major prefixes of the actual role matrices.
        widest = max(config["roles"], key=lambda row: row[1]*row[2])
        size = widest[1]*widest[2]
        require(size == config["flat_CRS_cells"], f"{name}: incorrect actual CRS size")
        group = []
        for i,cpu in enumerate(cpus,1):
            directory = args.runs/f"{name}-r{i}"
            row = read(directory/"result.json")
            require(row["passed"] and row["k"] == k and row["reference"] is False,
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
                       widest_actual_role=widest[0], widest_role_rows=widest[1], widest_role_columns=widest[2],
                       shared_crs_ring_elements=size, shared_crs_bytes=384*size,
                       backend_pivot_call_coefficient=prod(2*(1+v["s"]) for v in config["layers"]),
                       classical_log2=weak["log2_cost"]["core_classical"],
                       quantum_sensitivity_log2=weak["log2_cost"]["core_quantum"],
                       history_statistical_bits=-log2(history),
                       conditional_bits_if_cutoff_le_2m130_and_each_adv_le_2m132=-log2(conditional_budget))
        for key in ("prover_seconds","verifier_seconds","verify_replay_seconds","wall_seconds",
                    "prover_bytes","verifier_bytes","total_bytes","max_rss_kib"):
            values = [r[key] for r in group]
            for fn in (min,median,max):
                summary[key+"_"+fn.__name__] = fn(values)
        groups.append(summary)
    write_rows(out/"trials.csv",observations)
    write_rows(out/"summary.csv",groups)
    write_rows(out/"attacks.csv",attacks)
    (out/"wire-checks.json").write_text(json.dumps(wires,indent=2)+"\n",encoding="utf-8")
    lines = ["# Five-layer ordinary transitions",
             "",
             f"{machine['cpu']} ({machine['architecture']}), Linux `{machine['kernel']}`, "
             f"Python {machine['python']}. The three trials used logical CPUs {', '.join(map(str, cpus))}.",
             "",
             "| Profile | Prover (s) | Verifier (s) | Public replay (s) | Wall (s) | P→V (B) | V→P (B) | Total (B) | Peak RSS (MiB) |",
             "|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for r in groups:
        lines.append(f"| {r['profile']} | {r['prover_seconds_median']:.3f} | {r['verifier_seconds_median']:.3f} | "
                     f"{r['verify_replay_seconds_median']:.3f} | {r['wall_seconds_median']:.3f} | {r['prover_bytes_median']:,} | "
                     f"{r['verifier_bytes_median']:,} | {r['total_bytes_median']:,} | "
                     f"{r['max_rss_kib_max']/1024:.2f} |")
    lines += ["", "Times and wire lengths are medians of three fresh interactions. "
              "Peak RSS is the largest individual-process high-water mark across those trials. "
              "Common public inputs and CRS are available to both parties before timing the prover. "
              "Verifier time includes issuing and encoding coins plus witness-free verification of the decoded files. "
              "Tamper checking is timed separately. Wire lengths count the two actual padded streams, including retries. "
              "Wall time also covers fixture preparation, file I/O and the negative check. "
              "Each profile ran as a wave of three trials pinned to separate physical cores and L3 caches on a shared host.",
              "", "| Profile | Largest actual CRS role | Shared packed CRS (B) | Backend extractor pivot-call coefficient | Weakest generic role | Classical / quantum Core-SVP (log₂ work) | Issued-history statistical exponent |",
              "|---|---|---:|---:|---|---:|---:|"]
    for r in groups:
        lines.append(f"| {r['profile']} | {r['widest_actual_role']} ({r['widest_role_rows']} × {r['widest_role_columns']:,}) | "
                     f"{r['shared_crs_bytes']:,} | {r['backend_pivot_call_coefficient']:,} | {r['weakest_role']} | "
                     f"{r['classical_log2']:.3f} / {r['quantum_sensitivity_log2']:.3f} | {r['history_statistical_bits']:.4f} |")
    lines += [
              "", "The optimized profiles target 128-bit classical attack cost, with common screening floors of "
              "144 classical Core-SVP and 136 quantum Core-SVP sensitivity. The theorem assumes a classically rewindable prover "
              "and retains six time-bounded Module-SIS advantages and a global extraction cutoff. The issued-attempt cap is 16,384. "
              "The CRS column counts 384 encoded bytes per ring element. The extractor coefficient counts calls in a reduction.",
              "", "The [per-role costs](attacks.csv) label MATZOV extrapolations and identify finite CRT fits. "
              "Generic blocks above the fitted limit of 1,024 are " +
              ", ".join(f"{a['role']} at k={a['k']} (block {a['bkz_block']:,})" for a in attacks
                        if a['bkz_block'] > 1024) + "."]
    (out/"experiment-table.md").write_text("\n".join(lines)+"\n", encoding="utf-8")
    print(f"Verified {len(unique)} distinct verifier transcripts across {len(groups)} profiles.")


if __name__ == "__main__":
    main()
