"""Exhaustively check persistence of d-sortedness under an h-pass.

For ``2 <= n <= 8``, the program enumerates every permutation and every pair
of distances ``d`` and ``h``. Whenever the input is d-sorted, the output of a
standard gapped-insertion h-pass must remain d-sorted (``lem:persistence``).
"""
import itertools

from verification_utils import gapped_insertion, is_dsorted


def main():
    viol = []
    for n in range(2, 9):
        for d in range(1, n):
            for h in range(1, n):
                for perm in itertools.permutations(range(n)):
                    if not is_dsorted(perm, d):
                        continue
                    out, _ = gapped_insertion(perm, h)
                    if not is_dsorted(out, d):
                        viol.append((n, d, h, perm))
    print("Persistence violations:", len(viol))
    if viol:
        print(viol[:5])
    assert not viol
    print("PASS: exhaustive enumeration found no violations for n <= 8.")


if __name__ == "__main__":
    main()
