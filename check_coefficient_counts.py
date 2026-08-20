"""Check auxiliary coefficient-ball and positive-simplex counting formulas.

The program directly enumerates lattice points for ``1 <= s <= 4`` and
``0 <= t <= 6`` and compares them with the corresponding closed forms.
"""
import itertools
from math import comb


def V_s_t(s, t):
    """Return the closed form for the integer l1-ball cardinality."""
    return sum(2 ** j * comb(s, j) * comb(t, j) for j in range(min(s, t) + 1))


def brute_l1_ball(s, t):
    """Enumerate integer vectors with l1 norm at most ``t``."""
    cnt = 0
    for a in itertools.product(range(-t, t + 1), repeat=s):
        if sum(abs(x) for x in a) <= t:
            cnt += 1
    return cnt


def brute_pos_simplex(s, t):
    """Enumerate nonnegative integer vectors with coordinate sum at most ``t``."""
    cnt = 0
    for a in itertools.product(range(0, t + 1), repeat=s):
        if sum(a) <= t:
            cnt += 1
    return cnt


def main():
    bad1, bad2 = [], []
    for s in range(1, 5):
        for t in range(0, 7):
            if brute_l1_ball(s, t) != V_s_t(s, t):
                bad1.append((s, t, brute_l1_ball(s, t), V_s_t(s, t)))
            if brute_pos_simplex(s, t) != comb(s + t, s):
                bad2.append((s, t))
    print("Integer l1-ball mismatches:", bad1)
    print("Positive-simplex mismatches:", bad2)
    assert not bad1 and not bad2
    print("PASS: both counting formulas match direct enumeration.")


if __name__ == "__main__":
    main()
