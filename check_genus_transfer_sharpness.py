"""Check the explicit family making the genus-to-transfer bound sharp."""

from math import gcd

from verification_utils import apery_set, sg_reach


def exact_signed_length(A: int, B: int, h: int) -> int:
    best = None
    for u in range(-4 * h, 4 * h + 1):
        rem = h - u * A
        if rem % B:
            continue
        v = rem // B
        length = abs(u) + abs(v)
        best = length if best is None else min(best, length)
    assert best is not None
    return best


def main() -> int:
    checked = 0
    for h in range(3, 101):
        if h % 5 == 2:
            continue

        A, B = 2 * h + 1, 3 * h - 1
        assert gcd(A, B) == 1
        qgens = [5, A, B]
        conductor = 6 * h - 6
        M = 6 * h - 1
        qreach = sg_reach(qgens, M)
        gaps = [m for m in range(1, M + 1) if not qreach[m]]
        assert len(gaps) == 3 * h - 2
        assert max(gaps) == conductor - 1
        assert all(qreach[m] for m in range(conductor, M + 1))

        ap = apery_set(qgens, 5, M)
        assert sorted(ap.values()) == sorted([0, A, 2 * A, B, 2 * B])

        n = 6 * h * h
        prefix_reach = sg_reach([B, A], n)
        for m in range(M + 1):
            assert prefix_reach[m * h] == qreach[m]

        expected = h if h % 5 in (0, 4) else 2 * h
        assert exact_signed_length(A, B, h) == expected
        assert M >= 2 * len(gaps) + 1
        checked += 1

    print(f"Genus-transfer sharpness: parameters={checked}, violations=0")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
