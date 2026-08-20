"""Check extremal ray moments and finite-window ray-to-signed transfer."""

from math import gcd

from verification_utils import pos_len, random_gaps, ray_quotient


def main() -> int:
    moment_cases = transfer_cases = 0

    schedules = []
    for n in (24, 37, 53, 79):
        schedules.append((n, [5, 4, 1] if n == 24 else [n // 2, n // 3, 1]))
        schedules.extend(
            (n, random_gaps(n, seed=8100 + 100 * n + seed, lo=2, hi=7))
            for seed in range(30)
        )

    for n, H in schedules:
        H = sorted({h for h in H if 1 <= h < n}, reverse=True)
        if H[-1] != 1:
            H.append(1)
        for j in range(1, len(H) + 1):
            h = H[j - 1]
            M, Q, G, D, g = ray_quotient(H, n, j)
            rho = n - M * h

            lower = rho * g + h * g * (g - 1) // 2
            upper = n * g - h * g * (g + 1) // 2
            assert lower <= D <= upper
            lower_extremal = list(range(M - g + 1, M + 1)) if g else []
            upper_extremal = list(range(1, g + 1))
            assert (D == lower) == (G == lower_extremal)
            assert (D == upper) == (G == upper_extremal)
            # Exact integer form of g <= (1 + sqrt(1 + 8D/h))/2.
            assert h * g * (g - 1) <= 2 * D
            moment_cases += 1

            if j < 2 or M < 2 * g + 1:
                continue

            assert g >= 1 and 1 not in Q
            q = next(
                (m for m in range(1, 2 * g + 1) if m in Q and m + 1 in Q),
                None,
            )
            assert q is not None
            Lq = pos_len(q, H, j)
            Lq1 = pos_len(q + 1, H, j)
            assert Lq is not None and Lq1 is not None
            assert Lq <= q - 1 and Lq1 <= q
            assert Lq + Lq1 <= 4 * g - 1
            transfer_cases += 1

    # Regression for the counterexample that rules out the wrong ``-1`` root.
    n, H, j = 9, [5, 4, 1], 2
    _, _, _, D, g = ray_quotient(H, n, j)
    assert (g, D) == (2, 6)
    assert H[j - 1] * g * (g - 1) <= 2 * D

    print(
        "Ray moments and finite-window transfer: "
        f"moment_cases={moment_cases}, transfer_cases={transfer_cases}, violations=0"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
