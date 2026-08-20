# -*- coding: utf-8 -*-
"""
Independently check the transfer chains and numerical calibration.

Checks:
  1. Exact lambda_j (minimum signed l1 length) for the top six Pratt gaps at
     n = 10^6, plus the onset of the local lambda_j = 2 relation at the
     three scales used by the quantitative illustration.
  1b. Exact lambda_j = 2 for the Fibonacci-grid transfer chain used by
     the hybrid calibration (j = 3..12 at n = 10^6, 10^9, 10^12).
  2. Phase-separation search: does any character k satisfy
        ||k h_j / n||_T > 2L * delta_{j-1}(k)
     for the calibrated passes (Pratt h=54, Fibonacci h=13, n=1000, L=2)?
  3. Regenerate weighted rows (n=1000, Pratt/Fibonacci), supplementary profile
     summaries, and the n=36 worked-example values.
  4. Hybrid gain asymptotic: gain(n) = Theta(n^{1/3}) for fixed j0=2, s=6.
"""
import math
from fractions import Fraction


class SearchLimitExceeded(RuntimeError):
    """A finite signed representation exists, but exceeds the search limit."""


def torus_norm(x, n):
    """Exact torus norm ||x/n||_T."""
    r = x % n
    return Fraction(min(r, n - r), n)


def iroot(n, k):
    """Return max{x >= 0: x**k <= n} using exact integer arithmetic."""
    lo, hi = 0, 1
    while hi ** k <= n:
        hi *= 2
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if mid ** k <= n:
            lo = mid
        else:
            hi = mid
    return lo


def three_smooth_below(n):
    """All 2^a 3^b < n in decreasing order."""
    out = []
    a = 0
    while 2 ** a < n:
        b = 0
        while (2 ** a) * (3 ** b) < n:
            out.append((2 ** a) * (3 ** b))
            b += 1
        a += 1
    return sorted(set(out), reverse=True)


def fib_gaps(n=1000):
    f = [1, 1]
    while f[-1] < n:
        f.append(f[-1] + f[-2])
    g = [x for x in f if x < n]
    return sorted(set(g), reverse=True)


def gcd_list(xs):
    g = 0
    for x in xs:
        g = math.gcd(g, x)
    return g


def lam_exact(H, j, max_l1):
    """Return ``(lambda_j, witness)`` by exact bounded l1 minimization.

    ``lambda_j`` is ``math.inf`` exactly when the prefix gcd does not divide
    the target.  If a representation exists but its minimum exceeds
    ``max_l1``, raise ``SearchLimitExceeded`` instead of conflating that case
    with infinity.
    """
    target = H[j - 1]
    pref = H[: j - 1]
    g = gcd_list(pref)
    if target % g != 0:
        return math.inf, None

    # Divide out the common gcd before enumeration.  At a fixed l1 budget,
    # enumerate all but the final coefficient and solve the last one exactly.
    target //= g
    pref = [h // g for h in pref]
    d = len(pref)
    suffix_gcd = [0] * d
    suffix_max = [0] * d
    for i in range(d - 1, -1, -1):
        suffix_gcd[i] = math.gcd(pref[i], suffix_gcd[i + 1] if i + 1 < d else 0)
        suffix_max[i] = max(pref[i], suffix_max[i + 1] if i + 1 < d else 0)

    def search(pos, budget, residual, coeffs):
        if residual % suffix_gcd[pos] != 0:
            return None
        if abs(residual) > budget * suffix_max[pos]:
            return None
        if pos == d - 1:
            if residual % pref[pos] == 0:
                x = residual // pref[pos]
                if abs(x) <= budget:
                    return coeffs + [x]
            return None

        for x in range(-budget, budget + 1):
            witness = search(
                pos + 1,
                budget - abs(x),
                residual - x * pref[pos],
                coeffs + [x],
            )
            if witness is not None:
                return witness
        return None

    for length in range(max_l1 + 1):
        witness = search(0, length, target, [])
        if witness is not None:
            assert sum(x * h for x, h in zip(witness, H[: j - 1])) == H[j - 1]
            assert sum(abs(x) for x in witness) == length
            return length, witness

    raise SearchLimitExceeded(
        f"lambda_{j} is finite but exceeds the search limit {max_l1}"
    )


def check_top_pratt_lambdas():
    print("=== 1. exact lambda_j for the top Pratt prefix at n=10^6 ===")
    n = 10 ** 6
    H = three_smooth_below(n)
    expected = {2: math.inf, 3: math.inf, 4: 110, 5: math.inf, 6: 14}
    print(f"top gaps: {H[:6]}")
    for j, wanted in expected.items():
        lj, witness = lam_exact(H, j, max_l1=110)
        assert lj == wanted, f"lambda_{j}: expected {wanted}, got {lj}"
        shown = "inf" if math.isinf(lj) else str(lj)
        print(f"  j={j} h_j={H[j-1]}  lambda_j={shown}  witness={witness}")

    print("  onset of the local Pratt lambda_j=2 relation:")
    for n in (10 ** 6, 10 ** 9, 10 ** 12):
        H = three_smooth_below(n)
        j = next(idx for idx, h in enumerate(H, 1) if 3 * h < n)
        j0 = j - 1
        Nj0 = iroot(n - 1, j0)
        assert 3 * H[j - 1] < n
        assert j == 1 or 3 * H[j - 2] >= n
        assert Nj0 == 1
        print(
            f"    n=1e{int(math.log10(n))}: first transfer pass j={j}, "
            f"j0={j0}, N_j0={Nj0}"
        )


def check_fib_lambdas():
    print("\n=== 1b. exact lambda_j = 2 for the Fibonacci transfer chain ===")
    for n in (10 ** 6, 10 ** 9, 10 ** 12):
        H = fib_gaps(n)
        for j in range(3, 13):
            lj, witness = lam_exact(H, j, max_l1=2)
            assert lj == 2, f"n={n} j={j}: lambda_j={lj}, expected exactly 2"
            assert witness[j - 3] == 1 and witness[j - 2] == -1 and all(
                x == 0 for x in witness[: j - 3]
            ), f"n={n} j={j}: unexpected witness {witness}"
        print(
            f"  n=1e{int(math.log10(n))}: lambda_j = 2 exactly for j = 3..12 "
            "(witness h_j = h_{j-2} - h_{j-1})"
        )


def phase_separation_search():
    print("\n=== 2. phase-separation character search (n=1000, L=2) ===")
    n = 1000
    L = 2
    cases = [("Pratt h=54", three_smooth_below(n), 54),
             ("Fibonacci h=13", fib_gaps(n), 13)]
    for name, H, hj in cases:
        j = H.index(hj) + 1
        pref = H[: j - 1]
        best = None
        for k in range(1, n):
            d = max(torus_norm(k * h, n) for h in pref)
            t = torus_norm(k * hj, n)
            if t > 2 * L * d:
                best = (k, t, d, 2 * L * d)
                break
        if best:
            k, t, d, rhs = best
            print(
                f"  {name}: FOUND k={k}: ||kh_j/n||={float(t):.4f} "
                f"> 2L delta={float(rhs):.4f} (delta={float(d):.4f})"
            )
        else:
            print(
                f"  {name}: no k in 1..999 satisfies "
                f"||kh_j/n|| > 2L*delta_{{j-1}}(k)"
            )
        assert best is None, f"unexpected phase-separating character for {name}: {best}"


def shortest_positive_lengths(H, j, M):
    """Exact ell_j^->(m) for m in 0..M via positive-sum dynamic programming."""
    pref = H[: j - 1]
    target = H[j - 1]
    INF = float("inf")
    maxv = target * M
    best = [INF] * (maxv + 1)
    best[0] = 0
    for value in range(1, maxv + 1):
        best[value] = min(
            (best[value - h] + 1 for h in pref if h <= value),
            default=INF,
        )
    return [best[target * m] for m in range(M + 1)]


def check_weighted_profiles():
    print("\n=== 3. regenerate weighted rows and supplementary profiles (n=1000) ===")
    n = 1000
    Pratt = three_smooth_below(n)
    Fib = fib_gaps(n)
    expected = {
        ("Pratt", 162): (1, 838, 838),
        ("Pratt", 108): (1, 892, 892),
        ("Pratt", 81): (1, 919, 919),
        ("Pratt", 54): (1, 946, 946),
        ("Fibonacci", 610): (1, 390, 390),
        ("Fibonacci", 233): (4, 1670, 1670),
        ("Fibonacci", 89): (11, 5126, 5126),
        ("Fibonacci", 34): (29, 14210, 14210),
        ("Fibonacci", 13): (27, 22086, 20812),
    }
    actual = {}
    for name, H, hs in [("Pratt", Pratt, [162, 108, 81, 54]),
                        ("Fibonacci", Fib, [610, 233, 89, 34, 13])]:
        for h in hs:
            j = H.index(h) + 1
            M = (n - 1) // h
            pref = H[: j - 1]
            # full semigroup closure below n (iterate to fixpoint)
            reach = {0}
            changed = True
            while changed:
                changed = False
                for g in pref:
                    new = {s + g for s in reach if s + g < n} - reach
                    if new:
                        reach |= new
                        changed = True
            G = [m for m in range(1, M + 1) if m * h not in reach]
            D = sum(n - m * h for m in G)
            g = len(G)
            sharpened = n * g - h * g * (g + 1) // 2
            actual[(name, h)] = (g, sharpened, D)
            print(
                f"  {name:10s} h={h:4d}  g={g:3d}  "
                f"sharpened={sharpened:6d}  D={D}"
            )
    assert actual == expected, f"weighted-calibration mismatch: {actual}"

    print("  Pratt h=54 profile (m: ell):")
    j = Pratt.index(54) + 1
    prof = shortest_positive_lengths(Pratt, j, 18)
    print("   ", [(m, (int(x) if x != float("inf") else "inf")) for m, x in enumerate(prof)])
    assert [m for m, x in enumerate(prof) if x == float("inf")] == [1]
    assert prof[0] == 0
    assert max(x for x in prof if x != float("inf")) == 2

    print("  Fibonacci h=13 profile (m: ell), M=76:")
    j = Fib.index(13) + 1
    prof = shortest_positive_lengths(Fib, j, 76)
    missing = [m for m, x in enumerate(prof) if x == float("inf")]
    print(f"   missing count={len(missing)} missing={missing}")
    print("   (m:ell) =", [(m, int(x)) for m, x in enumerate(prof) if x != float("inf")])
    assert missing == list(range(1, 21)) + list(range(35, 42))
    assert prof[0] == 0 and prof[29] == 1 and prof[58] == 2 and prof[69] == 31
    assert max(x for x in prof if x != float("inf")) == 31

    print("  Worked example n=36, h_j=6:")
    H = [32, 27, 24, 18, 16, 12, 9, 8, 6]
    j = H.index(6) + 1
    reach = {0}
    changed = True
    while changed:
        changed = False
        for h in H[: j - 1]:
            new = {value + h for value in reach if value + h < 36} - reach
            if new:
                reach |= new
                changed = True
    Q = [m for m in range(6) if 6 * m in reach]
    G = [m for m in range(1, 6) if m not in Q]
    lj, witness = lam_exact(H, j, max_l1=2)
    assert 2 in Q and 3 in Q and G == [1]
    assert lj == 2 and witness is not None
    print(f"    Q contains 2,3; missing={G}; lambda_j={lj}; witness={witness}")


def gain_asymptotic():
    print("\n=== 4. hybrid gain for fixed j0=2, s=6 ===")
    expected_rounded = {10 ** 6: 6.9, 10 ** 9: 63.8, 10 ** 12: 631.3}
    for n in (10 ** 6, 10 ** 9, 10 ** 12):
        N2 = iroot(n - 1, 2)
        N6 = iroot(n - 1, 6)
        gain = Fraction(N2, (2 ** 4) * N6)
        assert round(float(gain), 1) == expected_rounded[n]
        print(
            f"  n=1e{int(math.log10(n))}: gain={float(gain):.1f}x  "
            f"(n^(1/3)/16 = {n ** (1/3) / 16:.3f})"
        )


def main():
    check_top_pratt_lambdas()
    check_fib_lambdas()
    phase_separation_search()
    check_weighted_profiles()
    gain_asymptotic()
    print("\nALL CALIBRATION CONSISTENCY CHECKS PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
