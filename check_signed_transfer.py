"""Check the signed-transfer inequality and its structural consequences.

The bounded breadth-first search is conservative: a window that is too small
can overestimate a finite transfer length, but cannot create a false short
representation.
"""
from verification_utils import (
    delta_max,
    fib_gaps,
    lam,
    local_apery_phase_bound,
    nrm,
    random_gaps,
    three_smooth,
)

def main():
    rng_seqs = []
    for n in (101, 128, 211, 256, 401):
        rng_seqs.append((three_smooth(n), n))
        rng_seqs.append((fib_gaps(n), n))
        rng_seqs.extend((random_gaps(n, seed=1000 + i), n) for i in range(6))

    tested = viol = finite = 0
    for H, n in rng_seqs:
        for j in range(2, len(H) + 1):
            lj, wit = lam(H, j)
            if lj is None:
                continue
            finite += 1
            assert sum(wit[i] * H[i] for i in range(j - 1)) == H[j - 1]
            assert sum(abs(t) for t in wit) == lj
            for k in range(1, min(n, 41)):
                tested += 1
                if nrm(k * H[j - 1], n) > lj * delta_max(H, j - 1, k, n):
                    viol += 1
    print(f"Signed-transfer inequality: tested={tested}, finite={finite}, violations={viol}")

    bad_ge2 = bad_gcd = 0
    from math import gcd
    for H, n in rng_seqs:
        for j in range(2, len(H) + 1):
            g = 0
            for i in range(j - 1):
                g = gcd(g, H[i])
            lj, _ = lam(H, j, window_mult=5)
            if lj is not None and lj < 2:
                bad_ge2 += 1
            if (H[j - 1] % g == 0) != (lj is not None):
                bad_gcd += 1
    inf2 = sum(1 for H, _ in rng_seqs if lam(H, 2)[0] is None)
    print(f"Lower bound lambda_j >= 2: violations={bad_ge2}")
    print(f"Gcd finiteness criterion: violations={bad_gcd}")
    print(f"Infinite lambda_2 cases: {inf2}/{len(rng_seqs)}")

    rows = []
    for H, n in [([21, 13, 8, 5, 3, 2, 1], 200),
                 ([34, 21, 13, 8, 5, 3, 2, 1], 300),
                 ([55, 34, 21, 13, 8, 5, 3, 2, 1], 400),
                 ([31, 17, 7, 1], 200), ([29, 23, 19, 1], 200)]:
        for j in range(3, len(H) + 1):
            lj, _ = lam(H, j)
            cb = local_apery_phase_bound(H, n, j)
            if lj is None or cb is None:
                continue
            rows.append((H, n, j, H[j - 1], lj, cb[0]))
    worse = sum(1 for r in rows if r[4] > r[5])
    best = max(rows, key=lambda r: r[5] / r[4])
    print(f"Transfer versus local Apéry bound: rows={len(rows)}, violations={worse}, "
          f"max_gain={best[5] / best[4]:.2f}x (H={best[0]}, j={best[2]})")

    return viol + bad_ge2 + bad_gcd + worse


if __name__ == "__main__":
    raise SystemExit(1 if main() else 0)
