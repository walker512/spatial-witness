# Spatial Witness

**Find missing reconstructed point samples, not just close-looking matches.**

A local, stdlib Python CLI that compares two **already aligned** point samples, records bidirectional nearest-point witnesses, and checks explicit global and regional coverage requirements. A reconstruction that kept only half a room can have perfect precision: the included example reports precision 2/2, recall 2/4, and a failed doorway region instead of calling that reconstruction complete.

This is the point-sample evaluation core of a spatial reconstruction quality checker. It **does not reconstruct a scene or validate physical truth**, and is not a full 3D viewer, automatic registration system, mesh inspector or AI model evaluator. Input can come from any authorized reconstruction workflow, AI-based or otherwise; the producer's identity is not verified.

## Run

Python 3.11+; no runtime dependencies. From this source directory:

```sh
python -m spatial_witness inspect examples/reference.csv --unit m
python -m spatial_witness compare examples/contract.json examples/reference.csv examples/complete.csv --out complete-review --require-pass
python -m spatial_witness verify complete-review --require-pass
python -m spatial_witness compare examples/contract.json examples/reference.csv examples/missing-doorway.csv --out missing-review --require-pass
# The last command writes the review and exits 3, as intended.
python -m spatial_witness verify missing-review
python -m unittest discover -s tests -v
```

Optional installation: `python -m pip install .`, then use `spatial-witness` in place of `python -m spatial_witness`. Packaging may obtain setuptools from your configured package source if not already installed. The application itself does not network, execute subprocesses, read environment secrets, or modify source files.

`compare` writes a **new** output directory containing the original `contract.json`, `reference.csv`, `reconstruction.csv`, and generated `report.json`, `distances.csv`, `review.md`. `verify` recomputes the three outputs from the bundled sources and requires byte equality. It ignores unrelated extra files; it does not certify the whole directory safe to share. A source or output change without a consistent rebuild fails verification. Someone able to replace all inputs and outputs can create a new consistent bundle: hashes are **not signatures, provenance, freshness or authenticity**.

Exit codes: 0 = comparison completed / bundle internally consistent (not necessarily quality passed); 2 = invalid input, changed output, unavailable output or other operational error; 3 = `--require-pass` and any declared requirement failed or could not be evaluated. `inspect` returns sample count, duplicate count and bounds only, with `quality_evaluated: false`. No reference means no quality comparison; empty reference is rejected.

## Input contract

CSV header must be exactly `x,y,z`, with three fixed-decimal coordinates per data row. Quoting follows normal CSV rules. No extra attributes, blank data rows, missing coordinates, NaN/Infinity, scientific notation or rounding. Use an explicit integer part (`0.5`, not `.5`). At most **1000 points per cloud**, 2 MiB per input file; this release targets selected patches and small samples, **not full resolution room scans**. Larger clouds are rejected, never silently subsampled. If you sample externally, record and justify your sampling policy separately: sample choice and density can change scores.

Each cloud declares `m` (up to 6 decimal places) or `mm` (up to 3 decimal places). Both normalize exactly to integer micrometres, bounded to ±1,000,000 metres per axis. A point's row index is its stable zero-based ID within that source. Duplicate rows are counted and **retained with full weight**, not silently deduplicated.

Contract JSON is demonstrated in `examples/contract.json`:

- `schema`: integer 1.
- `frame`: common declared coordinate-frame ID (ASCII letter followed by letters, digits, `_ . -`, up to 48 characters). **You** must correctly align both clouds and declare their units. This label does not prove alignment. There is no automatic rigid alignment, scale fitting, CRS conversion or ICP.
- `reference_unit`, `reconstruction_unit`: `m` or `mm`.
- `threshold_um`: positive integer micrometres, up to 10^12. A point matches only if **distance strictly less than threshold**. Equality is a miss.
- `minimum_precision_bp`, `minimum_recall_bp`: integer basis points, 0..10000. 9000 means 90%. Gates compare exact integer hit counts, not rounded display ratios.
- `regions`: 0..16 named inclusive axis-aligned boxes, with `min_um`, `max_um` (three bounded integer coordinates) and `minimum_recall_bp`. Regions can overlap. These test **reference samples inside each box against the entire reconstruction**, not reconstruction points cropped to that box. A sample just outside the box can match a reference point near its edge. Zero reference samples means `unknown_no_reference_samples`, never a pass, even if the minimum is zero.

JSON rejects unknown fields, duplicate keys, booleans in integer fields, floating values and nesting beyond 16 levels. Leaf input files must be regular and not symbolic links; special files are rejected before reading. This is a local trusted-directory tool, not a security boundary against hostile concurrent filesystem mutation or symlinked parent directories.

## Math and witnesses

For every reconstructed sample find the nearest reference sample; for every reference sample find the nearest reconstructed sample. Distances are squared Euclidean distances in integer µm². All comparisons are exact at the accepted input resolution. Equal-distance ties choose the first target row. The exhaustive search is O(N×M) per direction; no approximation or hidden registration is applied.

- Precision = matched reconstruction rows / all reconstruction rows.
- Recall = matched reference rows / all reference rows.
- F-score = harmonic mean of precision and recall, for presentation only.
- The report's ratios are floating-point presentation values. Counts and squared distances are the authoritative values; gates use exact integer arithmetic.
- Empty reconstruction: precision and F-score are undefined (`null`), recall is zero, nearest target indices/distances are `null`; overall status fails, even with zero gates.
- Overall `pass` requires defined global precision, global recall meeting their thresholds **and every region passing**. Otherwise status is `fail_or_unknown`; per-region status distinguishes unknown from failed.

`distances.csv` exposes direction, source row index, nearest target row, exact squared distance and strict match result. This helps locate missed or outlying samples in your own plotting/geometry tools. `report.json` includes the source-byte hashes, contract, summary bounds/duplicates, all witnesses and regional reference indices. All generated files use UTF-8/ASCII and LF bytes.

## Limits that matter

This compares **point samples**, not surfaces, topology, normals, textures, occluded geometry, camera poses, semantic correctness or mesh watertightness. Dense repeated points can dominate a score. There is no universal good threshold: choose tolerances and reference quality appropriate to your task. Synthetic examples prove the code's stated arithmetic, **not** real-world reconstruction quality, device accuracy or useful thresholds for your room. A bad reference or incorrect registration can produce misleading results. No ground truth means `inspect` only; the tool never invents one.

Reports are pre-serialized and bounded to 16 MiB per generated file before creating the output directory. Existing directories are never overwritten. A disk/write failure can leave an incomplete directory; it will not verify, and you should remove/review it before choosing a fresh output name. Verification is not an atomic transaction or signed receipt.

**Privacy:** bundled sources include the complete coordinates and region descriptions. Point clouds can reveal private room layouts. No anonymization, redaction, upload or automatic telemetry is performed. Review every file before sharing; use only authorized inputs. Error messages do not echo input values or local file paths.

## Existing tools and honest scope

[Open3D](https://www.open3d.org/docs/latest/tutorial/pipelines/icp_registration.html) provides registration and evaluation; [Tanks and Temples](https://www.tanksandtemples.org/tutorial/) defines established nearest-point precision, recall and F-score evaluation, with a [full evaluation toolbox](https://github.com/isl-org/TanksAndTemples/tree/master/python_toolbox/evaluation). [pointcloud_evaluation](https://github.com/pronif/pointcloud_evaluation) also evaluates reconstruction precision/recall. These offer capabilities beyond this bounded CLI, and the metrics here are **not novel**.

This project's practical integration is explicit unit/frame declarations, exact threshold decisions, regional sample-coverage gates and a portable recomputable witness bundle. This is not a claim those features are unique or that this tool beats existing workflows. Only primary documentation was inspected; **no competitor runtime, speed or accuracy comparison was performed**, and no official benchmark protocol conformance is claimed. No upstream code or dataset is bundled.

## Validation and license

33 public synthetic tests ran on Linux / Python 3.12, zero skipped: exact threshold boundary, mixed units, asymmetric precision/recall, empty inputs, duplicates, regional unknowns, basis-point gates, extreme coordinates, 100 seeded random cases against an independent straightforward squared-distance oracle, malformed input, FIFO rejection, LF bytes, changed bundles and CLI exit behavior. See `tests/test_spatial.py`. These are synthetic code checks, not measured real scans or user adoption. Windows/macOS and a new hosted CI run were not tested; symlink/FIFO tests have explicit platform skips where unsupported. No browser functionality is provided or claimed.

MIT, copyright walker512. All included fixtures are original synthetic coordinates. Please do not include secrets or private scans in public bug reports.
