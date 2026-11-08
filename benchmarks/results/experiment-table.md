# Five-layer ordinary transitions

| Profile | Prover (s) | Verifier (s) | Public replay (s) | Wall (s) | P→V (B) | V→P (B) | Total (B) | Peak RSS (MiB) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| k16 | 19.737 | 18.071 | 18.036 | 37.912 | 106,775 | 36,390 | 143,165 | 86.65 |
| k64 | 34.227 | 31.068 | 31.033 | 65.544 | 110,127 | 39,119 | 149,246 | 167.50 |
| k256 | 83.797 | 74.946 | 74.905 | 159.603 | 118,100 | 41,407 | 159,507 | 523.89 |
| k1024 | 258.324 | 226.150 | 226.093 | 487.709 | 126,145 | 44,216 | 170,361 | 2014.43 |

Times and wire lengths are medians of three fresh interactions. Peak RSS is the largest individual-process high-water mark across those trials. Common public inputs and CRS are available to both parties before timing the prover. Verifier time includes issuing and encoding coins plus witness-free verification of the decoded files. Tamper checking is timed separately. Wire lengths count the two actual padded streams, including retries. Wall time also covers fixture preparation, file I/O and the negative check. Each profile ran as a wave of three trials pinned to separate physical cores and L3 caches on a shared host.

| Profile | Largest actual CRS role | Shared packed CRS (B) | Backend extractor pivot-call coefficient | Weakest generic role | Classical / quantum Core-SVP (log₂ work) | Issued-history statistical exponent |
|---|---|---:|---:|---|---:|---:|
| k16 | A_3 (16 × 361) | 2,217,984 | 120,000 | A_1 | 156.220 / 141.775 | 131.3725 |
| k64 | A_0 (16 × 665) | 4,085,760 | 168,000 | A_0 | 150.672 / 136.740 | 131.2412 |
| k256 | A_0 (16 × 1,819) | 11,175,936 | 302,400 | A_0 | 171.696 / 155.820 | 131.0099 |
| k1024 | A_0 (16 × 4,733) | 29,079,552 | 786,240 | A_0 | 162.936 / 147.870 | 130.6355 |

The optimized profiles target 128-bit classical attack cost, with common screening floors of 144 classical Core-SVP and 136 quantum Core-SVP sensitivity. The theorem assumes a classically rewindable prover and retains six time-bounded Module-SIS advantages and a global extraction cutoff. The issued-attempt cap is 16,384. The CRS column counts 384 encoded bytes per ring element. The extractor coefficient counts calls in a reduction.

The [per-role costs](attacks.csv) label MATZOV extrapolations and identify finite CRT fits. Generic blocks above the fitted limit of 1,024 are A_aux at k=16 (block 1,359), A_aux at k=64 (block 1,131).
