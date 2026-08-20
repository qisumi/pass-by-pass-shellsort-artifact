"""Check standard genus and Frobenius identities for Apéry sets.

For several numerical semigroups, direct enumeration is compared with the
Apéry formulas used in ``thm:apery-bridge`` and ``prop:apery-exact``.
"""
from math import gcd

from verification_utils import apery_set, frobenius_of, genus_of


def main():
    families = [(3, 5), (4, 7), (5, 7, 11), (6, 9, 20), (7, 10),
                (11, 13, 17), (4, 6, 9)]
    bad1, bad2 = [], []
    for gens in families:
        if gcd(*gens) != 1:
            continue
        g = genus_of(gens)
        F = frobenius_of(gens)
        for a in gens:
            v = apery_set(gens, a)
            if any(x is None for x in v.values()):
                continue
            pred_g = sum(v.values()) / a - (a - 1) / 2
            pred_F = max(v.values()) - a
            if abs(pred_g - g) > 1e-9:
                bad1.append((gens, a, g, pred_g))
            if pred_F != F:
                bad2.append((gens, a, F, pred_F))
    print("Apéry genus mismatches:", bad1)
    print("Frobenius identity mismatches:", bad2)
    assert not bad1 and not bad2
    print("PASS: both Apéry identities match direct enumeration.")


if __name__ == "__main__":
    main()
