"""Check the rearrangement, phase-spectrum, and Fourier lower inequalities.

The randomized checks use fixed seeds. The Fourier certificate is also
checked by exhaustive enumeration of all two- and three-pass schedules for
``4 <= n <= 9``.
"""
import itertools
import math
import random
from math import pi, sqrt

from verification_utils import shellsort_counts, tnorm


def best_perm_value(x):
    """Return the maximum rearranged weighted sum."""
    n = len(x)
    idx = sorted(range(n), key=lambda i: x[i])
    lab = [0] * n
    for r, i in enumerate(idx):
        lab[i] = r + 1
    return sum(lab[i] * x[i] for i in range(n))


def check_rearrangement(trials=4000, seed=1):
    random.seed(seed)
    bad = 0
    for _ in range(trials):
        n = random.randint(2, 12)
        x = [random.uniform(-3, 3) for _ in range(n)]
        m = sum(x) / n
        x = [t - m for t in x]
        if best_perm_value(x) < (n / 4) * sum(abs(t) for t in x) - 1e-9:
            bad += 1
    print(f"Rearrangement check: trials={trials}, violations={bad}")
    return bad


def check_norm_comparison(trials=1500, seed=2):
    random.seed(seed)
    bad = []
    for _ in range(trials):
        n = random.randint(8, 200)
        s = random.randint(1, min(5, n - 1))
        H = random.sample(range(1, n), s)
        g = min(sum(math.sin(pi * k * h / n) ** 2 for h in H) / s
                for k in range(1, n))
        d2 = min(sqrt(sum(tnorm(k * h / n) ** 2 for h in H) / s)
                 for k in range(1, n))
        b = min(max(tnorm(k * h / n) for h in H) for k in range(1, n))
        if not (4 * d2 * d2 <= g + 1e-12 and g <= pi * pi * d2 * d2 + 1e-12):
            bad.append(("spectral bounds", n, H, g, d2))
        if not (d2 <= b + 1e-12 and b <= sqrt(s) * d2 + 1e-12):
            bad.append(("δ–β", n, H, b, d2))
        if not (b <= 0.5 * sqrt(s * g) + 1e-12):
            bad.append(("β ≤ ½√(sγ)", n, H, b, g))
    print(f"Phase-spectrum check: trials={trials}, violations={len(bad)}")
    for x in bad[:3]:
        print("  ", x)
    return len(bad)


def check_fourier_lb():
    bad = []
    for n in range(4, 10):
        for k in (2, 3):
            for c in itertools.combinations(range(1, n), k):
                H = sorted(c, reverse=True)
                if H[-1] != 1:
                    continue
                W = max(sum(shellsort_counts(p, H)[0])
                        for p in itertools.permutations(range(n)))
                for s in range(1, len(H)):
                    beta = min(max(tnorm(kk * H[j] / n) for j in range(s))
                               for kk in range(1, n))
                    rhs = n * n / (8 * pi * (H[s] + pi * n * beta))
                    if W < rhs - 1e-9:
                        bad.append((n, tuple(H), s, W, rhs, beta))
    print("Fourier-certificate exhaustive violations:", len(bad))
    for x in bad[:5]:
        print("  ", x)
    return len(bad)


def two_pass_quadratic():
    print("Exact worst-case exchanges for two-pass schedules (h, 1):")
    for n in range(6, 9):
        row = []
        for h in range(2, n):
            W = max(sum(shellsort_counts(p, [h, 1])[0])
                    for p in itertools.permutations(range(n)))
            row.append(f"h={h}:{W}")
        print(f"  n={n}: " + "  ".join(row))


def main():
    bad = 0
    bad += check_rearrangement()
    bad += check_norm_comparison()
    bad += check_fourier_lb()
    two_pass_quadratic()
    assert bad == 0
    print("PASS: all Fourier-certificate checks succeeded.")


if __name__ == "__main__":
    main()
