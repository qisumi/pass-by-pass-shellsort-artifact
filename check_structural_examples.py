"""Reproduce the coefficient-entropy and length-genus counterexamples."""
from verification_utils import lam, sg_reach


def check_entropy_nogo():
    print("Coefficient-entropy example: conductor and genus")
    for M in (8, 12, 16, 20):
        r = sg_reach([M, M + 1, 2 * M + 1], 4 * M * M)
        F = max(x for x in range(1, 4 * M * M) if not r[x])
        g = sum(1 for x in range(1, 4 * M * M) if not r[x])
        ok = (F + 1 == M * (M - 1)) and (g == M * (M - 1) // 2)
        print(f"  M={M}: conductor={F+1} (expected {M*(M-1)}), genus={g} "
              f"(expected {M*(M-1)//2}) {'OK' if ok else 'MISMATCH'}")
        assert ok


def check_length_genus_nogo():
    print("Length-genus example: H=(q+1, q, 1), n=2q")
    for q in (3, 5, 8, 11):
        n = 2 * q
        H = [q + 1, q, 1]
        S2 = sg_reach(H[:2], n)
        M3 = n - 1
        Q3pos = [m for m in range(1, M3 + 1) if S2[m]]
        g3 = M3 - len(Q3pos)
        lambda3, witness = lam(H, 3)
        ok = (
            Q3pos == [q, q + 1]
            and g3 == 2 * q - 3
            and lambda3 == 2
            and sum(witness[i] * H[i] for i in range(2)) == 1
        )
        print(f"  q={q}: Q_3∩[1,M_3]={Q3pos} (expected [{q},{q+1}]), "
              f"g_3={g3} (expected {2*q-3}), lambda_3={lambda3} "
              f"(expected 2) {'OK' if ok else 'MISMATCH'}")
        assert ok


def main():
    check_entropy_nogo()
    check_length_genus_nogo()
    print("PASS: both structural examples match exactly.")


if __name__ == "__main__":
    main()
