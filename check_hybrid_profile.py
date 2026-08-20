"""Check signed propagation and the hybrid approximation-transfer profile.

The program checks the pointwise transfer inequality, its prefix consequence,
the chained hybrid bound, and compatibility with the classical Pratt upper
scale. Exact integer roots and rational arithmetic are used where possible.
"""
from fractions import Fraction
import math

from verification_utils import (
    beta_exact,
    delta_max,
    fib_gaps,
    lam,
    nrm,
    random_gaps,
    three_smooth,
)


def iroot(n, k):
    """Return the largest integer ``x`` such that ``x**k <= n``."""
    lo, hi = 0, 1
    while hi ** k <= n:
        hi *= 2
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if mid ** k <= n:
            lo = mid
        else:
            hi = mid
    return lo


def main():
    cases = []
    for n in (101, 151, 211, 257):
        cases.append((three_smooth(n), n))
        cases.append((fib_gaps(n), n))
        cases.extend((random_gaps(n, seed=2000 + i, lo=3, hi=6), n) for i in range(5))

    tested = viol = 0
    for H, n in cases:
        for s in range(1, len(H)):
            lj, _ = lam(H, s + 1)
            if lj is None:
                continue
            bs, _ = beta_exact(H, n, s)
            bs1, _ = beta_exact(H, n, s + 1)
            tested += 1
            if bs1 > lj * bs:
                viol += 1
    print(f"Prefix propagation: tested={tested}, violations={viol}")

    tested2 = viol2 = 0
    for H, n in cases:
        for s in range(1, len(H)):
            lj, _ = lam(H, s + 1)
            if lj is None:
                continue
            for k in range(1, min(n, 31)):
                tested2 += 1
                if nrm(k * H[s], n) > lj * delta_max(H, s, k, n):
                    viol2 += 1
    print(f"Pointwise propagation: tested={tested2}, violations={viol2}")

    tot = bad = 0
    for n in (201, 401, 601):
        H = three_smooth(n)
        for j in range(2, len(H) + 1):
            if 3 * H[j - 1] < n:
                tot += 1
                if lam(H, j)[0] != 2:
                    bad += 1
    print(f"Pratt lambda_j=2 relation: gaps={tot}, mismatches={bad}")

    tested3 = viol3 = 0
    for H, n in cases:
        for s in range(2, len(H)):
            bs, _ = beta_exact(H, n, s)
            for j0 in range(2, s + 1):
                bj0, _ = beta_exact(H, n, j0)
                prod = Fraction(1)
                ok = True
                for j in range(j0 + 1, s + 1):
                    lj, _ = lam(H, j)
                    if lj is None:
                        ok = False
                        break
                    prod *= lj
                if not ok:
                    continue
                tested3 += 1
                if bs > bj0 * prod:
                    viol3 += 1
    print(f"Hybrid chain: tested={tested3}, violations={viol3}")

    print("Hybrid-profile gain for the exact Fibonacci lambda_j=2 chain:")
    print("        n      s   Dirichlet     HYBRID       gain   j0")
    for n in (10 ** 6, 10 ** 9, 10 ** 12):
        for s in (6, 12):
            Ns = iroot(n - 1, s)
            dirich = Fraction(1, Ns)
            best = min(((Fraction(2 ** (s - j0), iroot(n - 1, j0)), j0)
                        for j0 in range(2, s + 1)
                        if iroot(n - 1, j0) >= 1),
                       key=lambda t: min(Fraction(1, 2), t[0]))
            B, j0 = best
            B = min(Fraction(1, 2), B)
            print(f"      1e{int(math.log10(n)):<3d} {s:4d}   {float(dirich):.3e}   "
                  f"{float(B):.3e}   {float(dirich / B):7.1f}x  j0={j0}")

    print("Compatibility check: hybrid lower bound / (n log^2 n)")
    bad_cons = 0
    for n in (10 ** 4, 10 ** 6, 10 ** 8, 10 ** 10, 10 ** 12):
        H = three_smooth(n)
        p = len(H)
        best = 0.0
        for s in range(2, p):
            hs1 = H[s] if s < p else 1
            for j0 in range(2, s + 1):
                Nj = iroot(n - 1, j0)
                if Nj < 1:
                    continue
                beta = Fraction(2 ** (s - j0), Nj)
                if beta > Fraction(1, 2):
                    continue
                best = max(best, n * n / (8 * math.pi * (hs1 + math.pi * n * float(beta))))
        ratio = best / (n * math.log2(n) ** 2)
        if ratio >= 1:
            bad_cons += 1
        print(f"      n=1e{int(math.log10(n)):<3d} ratio={ratio:.2e}  "
              f"{'OK' if ratio < 1 else 'CONTRADICTION'}")

    return viol + viol2 + bad + viol3 + bad_cons


if __name__ == "__main__":
    raise SystemExit(1 if main() else 0)
