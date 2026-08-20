"""Check the weighted ray-defect refinement and exact Apéry evaluation.

The checks cover the truncated-genus refinement, the Selmer and Brown-Shiue
formulas, and the saturated identity ``D_j = n*g(Q_j) - h_j*n_1(Q_j)``.
"""
from math import gcd

from verification_utils import (
    apery_genus_sylvester,
    fib_gaps,
    random_gaps,
    ray_quotient,
    sg_reach,
    three_smooth,
)


def main():
    tested = viol = neg = 0
    gains = []
    for n in (200, 401, 601, 1001):
        pool = [three_smooth(n), fib_gaps(n)]
        pool += [random_gaps(n, seed=3000 + i, lo=2, hi=8) for i in range(20)]
        for H in pool:
            for j in range(1, len(H) + 1):
                M, Q, G, D, g = ray_quotient(H, n, j)
                if g == 0:
                    continue
                tested += 1
                assert D == n * g - H[j - 1] * sum(G), "weighted-defect identity failed"
                if sum(G) < g * (g + 1) // 2:
                    viol += 1
                new = n * g - H[j - 1] * g * (g + 1) // 2
                if new < 0:
                    neg += 1
                if D > new:
                    viol += 1
                if new > 0:
                    gains.append((n * g) / new)
    gs = sorted(gains)
    print(f"Truncated-genus refinement: tested={tested}, violations={viol}, "
          f"negative_RHS={neg} median_gain={gs[len(gs) // 2]:.2f}x max={max(gs):.2f}x")

    bad_g = bad_n1 = chk = 0
    for gens in [(2, 3), (3, 5), (3, 7), (4, 5), (5, 7), (4, 6, 7),
                 (5, 6, 9), (6, 7, 8), (7, 11, 13), (3, 11), (2, 9), (10, 11)]:
        gg = 0
        for x in gens:
            gg = gcd(gg, x)
        if gg != 1:
            continue
        N = 400
        res = apery_genus_sylvester(gens, gens[0], N)
        if res is None:
            continue
        g_f, n1_f = res
        reach = sg_reach(gens, N)
        gaps = [x for x in range(1, N + 1) if not reach[x]]
        chk += 1
        if g_f != len(gaps):
            bad_g += 1
        if n1_f != sum(gaps):
            bad_n1 += 1
    print(f"Apéry moment formulas: tested={chk}, genus_mismatches={bad_g}, "
          f"first-moment mismatches={bad_n1}")

    appl = bad = 0
    for n in (200, 401, 601, 1001):
        pool = [three_smooth(n)]
        pool += [random_gaps(n, seed=4000 + i, lo=2, hi=7) for i in range(20)]
        for H in pool:
            for j in range(2, len(H) + 1):
                M, Q, G, D, g = ray_quotient(H, n, j)
                Qpos = [m for m in Q if m > 0]
                if not Qpos:
                    continue
                gg = 0
                for x in Qpos:
                    gg = gcd(gg, x)
                if gg != 1:
                    continue
                reach = sg_reach(Qpos, M)
                gapsQ = [m for m in range(1, M + 1) if not reach[m]]
                # The truncated window must contain every gap of Q_j.
                if gapsQ and max(gapsQ) > M:
                    continue
                appl += 1
                if D != n * len(gapsQ) - H[j - 1] * sum(gapsQ):
                    bad += 1
    print(f"Saturated weighted identity: applicable={appl}, mismatches={bad}")

    return viol + bad_g + bad_n1 + bad


if __name__ == "__main__":
    raise SystemExit(1 if main() else 0)
