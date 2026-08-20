"""Run all checks for the pass-by-pass Shellsort artifact."""

from pathlib import Path
import subprocess
import sys


CHECKS = (
    "check_coefficient_counts.py",
    "check_persistence.py",
    "check_apery_identities.py",
    "check_directed_upper_certificate.py",
    "check_fourier_certificate.py",
    "check_local_apery_transference.py",
    "check_structural_examples.py",
    "check_signed_transfer.py",
    "check_hybrid_profile.py",
    "check_weighted_apery_evaluation.py",
    "check_ray_signed_transference.py",
    "check_affine_ray_separation.py",
    "check_genus_transfer_sharpness.py",
    "reproduce_quantitative_calibration.py",
    "check_calibration_consistency.py",
)


def main() -> int:
    root = Path(__file__).resolve().parent
    failures = []

    for filename in CHECKS:
        print("=" * 72, flush=True)
        print(f"Running {filename}", flush=True)
        result = subprocess.run([sys.executable, str(root / filename)], check=False)
        if result.returncode:
            failures.append((filename, result.returncode))

    print("=" * 72)
    if failures:
        print("Verification failed:")
        for filename, returncode in failures:
            print(f"  {filename}: exit code {returncode}")
        return 1

    print(f"All {len(CHECKS)} verification programs passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
