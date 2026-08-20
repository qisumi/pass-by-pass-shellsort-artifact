"""Check the affine global-to-ray separation family exactly."""

from verification_utils import sg_reach


E = {1, 2, 3, 6, 7, 11}


def formula_member(x: int, t: int) -> bool:
    q, r = divmod(x, t)
    y = q - 6 * r
    if y < 0:
        return False
    return y not in E


def main() -> int:
    checked = 0
    for t in range(3, 101):
        n = 16 * t * t
        conductor = 6 * t * t + 6 * t
        gaps = [x for x in range(1, conductor) if not formula_member(x, t)]
        assert len(gaps) == 3 * t * t + 3 * t
        assert max(gaps) == conductor - 1
        assert sum(gaps) == t * (12 * t**3 + 19 * t**2 + 33 * t - 4) // 2
        assert all(formula_member(x, t) for x in range(conductor, conductor + 4 * t))

        if t <= 30:
            reach = sg_reach([4 * t, 5 * t, 6 * t + 1], conductor + 4 * t)
            for x in range(conductor + 4 * t + 1):
                assert reach[x] == formula_member(x, t)

        M = (n - 1) // (3 * t)
        ray_gaps = [m for m in range(1, M + 1) if not formula_member(3 * t * m, t)]
        assert ray_gaps == [1, 2]
        D = sum(n - 3 * t * m for m in ray_gaps)
        assert D == 32 * t * t - 9 * t

        earlier = [4 * t, 5 * t, 6 * t + 1]
        signed_summands = earlier + [-x for x in earlier]
        assert 3 * t not in signed_summands
        assert all(a + b != 3 * t for a in signed_summands for b in signed_summands)
        assert 2 * (4 * t) - 5 * t == 3 * t

        U = n * len(gaps) - sum(gaps)
        assert U == t * (84 * t**3 + 77 * t**2 - 33 * t + 4) // 2
        checked += 1

    print(f"Affine ray separation: parameters={checked}, violations=0")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
