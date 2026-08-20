# Computational checks for pass-by-pass Shellsort bounds

This repository contains the computational artifact accompanying the paper
*Which Gap Relations Matter, and When? Pass-by-Pass Upper and Lower Bounds for
Shellsort*.

Repository: <https://github.com/qisumi/pass-by-pass-shellsort-artifact>

The programs provide finite exhaustive checks and deterministic representative
family calculations for the paper's order-sensitive upper and lower bounds.
They are consistency checks, not computer-assisted proofs: every theorem is
proved analytically in the paper.

## Requirements

- Python 3.9 or newer
- No third-party packages

Randomized instance generators use fixed seeds.

## Quick start

From the repository root, run the complete suite:

```bash
python run_verification.py
```

Each check can also be run independently, for example:

```bash
python check_directed_upper_certificate.py
python check_fourier_certificate.py
```

The complete suite may take several minutes because some programs exhaustively
enumerate small Shellsort instances.

## Contents

The manuscript references use stable LaTeX labels rather than theorem numbers,
which may change during typesetting.

| Program | Manuscript labels | Purpose |
|---|---|---|
| `check_coefficient_counts.py` | -- | Checks auxiliary signed-ball and nonnegative-simplex counting identities. |
| `check_persistence.py` | `lem:persistence` | Exhaustively checks persistence of already-sorted distances. |
| `check_apery_identities.py` | `thm:apery-bridge`, `prop:apery-exact` | Checks the Apéry genus and Frobenius identities. |
| `check_directed_upper_certificate.py` | `thm:directed-upper` | Exhaustively checks the passwise and prefix-cutoff upper bounds. |
| `check_fourier_certificate.py` | `lem:rearrangement`, `thm:fourier-lb` | Checks the rearrangement, phase-spectrum, and Fourier lower inequalities. |
| `check_local_apery_transference.py` | `thm:apery-bridge` | Checks conductor and genus consequences of local Apéry witnesses. |
| `check_structural_examples.py` | `prop:length-genus-nogo` | Reproduces the length-genus counterexample and an auxiliary entropy example. |
| `check_signed_transfer.py` | `lem:signed-transfer`, `thm:apery-bridge` | Checks signed-transfer inequalities and exact transfer lengths. |
| `check_hybrid_profile.py` | `prop:hybrid-profile` | Checks propagation and hybrid approximation profiles. |
| `check_weighted_apery_evaluation.py` | `prop:truncated-genus-refinement`, `prop:apery-exact` | Checks weighted ray defects and exact Apéry evaluation. |
| `check_ray_signed_transference.py` | `prop:truncated-genus-refinement`, `thm:ray-genus-transfer`, `cor:defect-transfer` | Checks extremal moment bounds and finite-window genus-to-transfer constructions. |
| `check_affine_ray_separation.py` | `prop:affine-ray-separation` | Checks the affine global-to-pass separation family. |
| `check_genus_transfer_sharpness.py` | `prop:genus-transfer-sharpness` | Checks the explicit sharpness family. |
| `reproduce_quantitative_calibration.py` | `lem:dirichlet`, `prop:hybrid-profile`, `prop:truncated-genus-refinement`, `app:verification` | Reproduces the deterministic numerical calibration. |
| `check_calibration_consistency.py` | `lem:signed-transfer`, `prop:hybrid-profile`, `thm:short-ray-phase`, `prop:truncated-genus-refinement`, `app:verification` | Independently checks transfer chains, phase separation, weighted profiles, and worked examples. |
| `verification_utils.py` | -- | Shared exact-arithmetic, numerical-semigroup, and Shellsort utilities. |

## Reproducibility scope

The suite exhausts persistence for `2 <= n <= 8`, the directed upper bound for
decreasing two- and three-pass schedules ending in `1` with `4 <= n <= 8`, and
the Fourier inequalities for those schedules with `4 <= n <= 9`. The remaining
programs use deterministic families of signed transfers, local Apéry witnesses,
hybrid propagations, saturated Apéry evaluations, and weighted ray defects.
Individual programs print their tested case counts and other scope details.

The fixed seeds are recorded in the relevant source files. In particular, the
Fourier checks use seeds `1` and `2`; the local Apéry check uses seed `7`;
signed-transfer schedules use seeds `1000`--`1005`; hybrid schedules use seeds
`2000`--`2004`; weighted checks use seeds `3000`--`3019` and `4000`--`4019`;
and ray-transfer schedules use `8100 + 100*n + seed`.

The asymptotic results are not inferred from these finite computations.

The companion paper has a separate artifact at
<https://github.com/qisumi/full-product-shellsort-artifact>. The persistence
check is intentionally included in both repositories because both papers state
and use the same elementary lemma. The remaining programs are separated by
paper.

## Repository layout

All programs are kept at the repository root so they can be inspected and run
without package installation.

## License

No license is included in this revision. The source is publicly available for
inspection and reproducibility, but reuse rights remain reserved unless a
license is added later.
