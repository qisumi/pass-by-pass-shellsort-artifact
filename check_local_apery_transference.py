"""Check the local Apéry transference bounds on seeded semigroup instances.

For each numerical semigroup, shortest nonnegative representations determine
``L`` and ``B``. Directly computed conductors and genera are compared with
``cond(Q) <= L*B`` and ``g(Q) < L*B`` from ``thm:apery-bridge``.
"""
import collections
import random
from math import gcd

from verification_utils import sg_reach


def shortest_lengths(gens, N):
    """Return shortest nonnegative representation lengths up to ``N``."""
    INF = float("inf")
    L1 = [INF] * (N + 1)
    L1[0] = 0
    dq = collections.deque([0])
    while dq:
        x = dq.popleft()
        for g in gens:
            if x + g <= N and L1[x + g] == INF:
                L1[x + g] = L1[x] + 1
                dq.append(x + g)
    return L1


def witnesses(gens, a, N):
    """Return the least reachable value in every residue class modulo ``a``."""
    r = sg_reach(gens, N)
    ws = []
    for res in range(a):
        m = None
        x = res
        while x <= N:
            if r[x]:
                m = x
                break
            x += a
        if m is None:
            return None
        ws.append(m)
    return ws


def main():
    random.seed(7)
    N = 8000
    bad = tested = 0
    for _ in range(4000):
        k = random.randint(2, 4)
        gens = sorted(random.sample(range(2, 40), k))
        if gcd(*gens) != 1:
            continue
        r = sg_reach(gens, N)
        F = max(x for x in range(1, N) if not r[x])
        cond = F + 1
        genus = sum(1 for x in range(1, N) if not r[x])
        L1 = shortest_lengths(gens, N)
        B = max(gens)
        for a in gens:
            if a < 2:
                continue
            ws = witnesses(gens, a, N)
            if ws is None:
                continue
            L = max(L1[w] for w in ws)
            L = max(L, L1[a])
            LB = L * B
            tested += 1
            if not (cond <= LB and genus < LB):
                bad += 1
                if bad <= 3:
                    print("  mismatch", gens, "a=", a, "L=", L, "B=", B,
                          "LB=", LB, "cond=", cond, "genus=", genus)
    print(f"Local Apéry instances: tested={tested}, violations={bad}")
    assert bad == 0
    print("PASS: conductor and genus bounds hold on every seeded instance.")


if __name__ == "__main__":
    main()
