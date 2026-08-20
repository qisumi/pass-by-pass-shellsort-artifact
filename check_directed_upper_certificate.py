"""Exhaustively check the directed-semigroup upper certificate.

For every decreasing two- or three-pass schedule ending in 1, every
permutation is tested for ``4 <= n <= 8``. The program checks the passwise
bound ``T_j <= D_j`` and the prefix-cutoff bound defining ``mathfrak U_n`` in
``thm:directed-upper``.
"""
import itertools

from verification_utils import D_j, U_s, shellsort_counts


def main():
    viol = []
    slack = []
    seqs = []
    for n in range(4, 9):
        for k in (2, 3):
            for c in itertools.combinations(range(1, n), k):
                H = sorted(c, reverse=True)
                if H[-1] != 1:
                    continue
                seqs.append((n, tuple(H)))

    for n, H in seqs:
        D = [D_j(n, H, j) for j in range(1, len(H) + 1)]
        worst = 0
        for perm in itertools.permutations(range(n)):
            Ts, out = shellsort_counts(perm, H)
            assert out == sorted(out), "Shellsort implementation did not sort the input"
            for j, d in enumerate(D):
                if Ts[j] > d:
                    viol.append(("pass", n, H, j + 1, Ts[j], d))
            worst = max(worst, sum(Ts))
        U = min(sum(D[:s]) + U_s(n, H, s) for s in range(len(H) + 1))
        if worst > U:
            viol.append(("cutoff", n, H, worst, U))
        slack.append((n, H, worst, U, sum(D)))

    print("Directed-certificate violations:", len(viol))
    for v in viol[:8]:
        print("  ", v)
    print("Samples (n, H, exact worst exchanges, cutoff bound, sum D_j):")
    for row in slack[:8]:
        print("  ", row)
    assert not viol
    print("PASS: exhaustive enumeration found no violations for n <= 8.")


if __name__ == "__main__":
    main()
