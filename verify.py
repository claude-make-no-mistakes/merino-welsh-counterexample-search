"""Check native values by deletion–contraction and the matrix-tree theorem."""

import argparse
from collections import Counter
from fractions import Fraction
from functools import cache
from itertools import combinations
import json
from pathlib import Path
import random
import re
import subprocess


def check(condition, failure):
    if not condition:
        raise ValueError(failure)


def expand(n, bundles, mask):
    edges = []
    for i, (u, v, copies) in enumerate(bundles):
        doubled = 0 if mask & (1 << i) else 1
        if doubled:
            edges.extend([(u, v), (u, v)])
        for _ in range(copies - doubled):
            edges.extend([(u, n), (v, n)])
            n += 1
    return n, tuple(sorted(edges))


def determinant(matrix):
    previous = 1
    for i in range(len(matrix)-1):
        pivot = matrix[i][i]
        check(pivot > 0, "Laplacian cofactor must be positive definite")
        for j in range(i+1, len(matrix)):
            for k in range(i+1, len(matrix)):
                numerator = pivot * matrix[j][k] - matrix[j][i] * matrix[i][k]
                value, remainder = divmod(numerator, previous)
                check(remainder == 0, "Bareiss division was not exact")
                matrix[j][k] = value
        previous = pivot
    return matrix[-1][-1]


def trees(n, edges):
    minor = [[0]*(n-1) for _ in range(n-1)]
    for u, v in edges:
        minor[u][u] += 1
        if v != n-1:
            minor[v][v] += 1
            minor[u][v] -= 1
            minor[v][u] -= 1
    return determinant(minor)


def tutte(edges):
    """Direct deletion–contraction of the expanded graph, with no bundle identities."""
    @cache
    def visit(graph):
        if not graph:
            return 1, 1, 1
        for i, (u, v) in enumerate(graph):
            if u == v:
                _, c, t = visit(graph[:i] + graph[i+1:])
                return 0, 2*c, t

        degree = Counter(v for edge in graph for v in edge)
        copies = Counter(graph)
        e = min(copies, key=lambda edge: (-copies[edge], min(degree[v] for v in edge)))
        u, v = e
        i = graph.index(e)
        deletion = graph[:i] + graph[i+1:]
        neighbors = {w: set() for w in degree}
        for a, b in deletion:
            neighbors[a].add(b)
            neighbors[b].add(a)
        reached, stack = {u}, [u]
        while stack:
            for w in neighbors[stack.pop()] - reached:
                stack.append(w)
                reached.add(w)
        contracted = {w: w - (w > v) for w in degree}
        contracted[v] = u
        contraction = tuple(sorted(tuple(sorted((contracted[a], contracted[b])))
                                   for a, b in deletion))
        a_c, c_c, t_c = visit(contraction)
        if v not in reached:
            return 2*a_c, 0, t_c
        a_d, c_d, t_d = visit(deletion)
        return a_c+a_d, c_c+c_d, t_c+t_d

    result = visit(edges)
    visit.cache_clear()
    return result


def controls(search):
    rng = random.Random(6261125)
    examples = [(2, [(0, 1, k)]) for k in range(2, 8)]
    examples.append((4, [(u, v, 1) for u, v in combinations(range(4), 2)]))
    for _ in range(30):
        n = rng.randint(3, 5)
        support = {(u, u+1) for u in range(n-1)} | {(0, n-1)}
        support |= set(rng.sample(list(combinations(range(n), 2)), 2))
        bundles = [(u, v, 1) for u, v in sorted(support)]
        for i in rng.sample(range(len(bundles)), min(2, 9-len(bundles))):
            u, v, _ = bundles[i]
            bundles[i] = u, v, rng.randint(2, 3)
        examples.append((n, bundles))
    verified = 0
    for n, bundles in examples:
        header = [n, len(bundles), *(number for bundle in bundles for number in bundle)]
        text = " ".join(map(str, header)) + "\n"
        result = subprocess.run([search, "--table"], input=text.encode(),
                                capture_output=True, check=True)
        table = [json.loads(line) for line in result.stdout.splitlines()]
        check(len(table) == 1 << len(bundles), "Incomplete native table")
        for index, (mask, *expected) in enumerate(table):
            expanded_n, edges = expand(n, bundles, mask)
            check(index == mask and tutte(edges) == tuple(expected), "Tutte control differs")
            check(trees(expanded_n, edges) == expected[2], "Tree control differs")
            verified += 1
    return verified


def generator_count(n, m, geng, multig):
    if n == 2:
        return 1
    with subprocess.Popen([geng, "-q", "-C", str(n), f"{n}:{m}"],
                          stdout=subprocess.PIPE) as simple:
        counted = subprocess.run([multig, "-u", f"-e{m}:{m}"], stdin=simple.stdout,
                                 capture_output=True, check=True)
        simple.stdout.close()
        check(simple.wait() == 0, "Support generator failed")
    match = re.search(rb"(\d+) multigraphs generated", counted.stderr)
    check(match is not None, "Missing generator count")
    return int(match[1])


def verify(ledgers, search, generators):
    seen = set()
    totals = Counter()
    best = None
    worst = None
    for ledger in ledgers:
        for line in ledger.read_text(encoding="utf-8").splitlines():
            row = json.loads(line)
            n, m = row["vertices"], row["core_edges"]
            check((n, m) not in seen, "Duplicate stratum")
            seen.add((n, m))
            if generators is not None:
                check(generator_count(n, m, *generators) == row["cores"],
                      f"Generator count differs at order {n}, size {m}")
            totals["cores"] += row["cores"]
            totals["series_parallel_cores"] += row["series_parallel_cores"]
            totals["decorations"] += row["decorations"]
            totals["failures"] += row["failures"]
            case = row["minimum"]
            if case is None:
                check(row["series_parallel_cores"] == row["cores"] and not row["decorations"],
                      "Missing nonspecial cases")
                continue
            bundles, mask = case["bundles"], case["series_mask"]
            check(sum(edge[2] for edge in bundles) == m, "Wrong core size")
            expanded_n, edges = expand(n, bundles, mask)
            check(trees(expanded_n, edges) == case["trees"], "Tree count failed")
            ratio = Fraction(case["acyclic"] * case["cyclic"], case["trees"]**2)
            record = (ratio, edges, case)
            if best is None or ratio < best[0]:
                best = record
            if worst is None or m > worst[0]:
                worst = (m, record)
    check(best is not None, "Expected at least one non-series-parallel core")
    bound = max(m for _, m in seen)
    check(seen == {(n, m) for m in range(2, bound+1) for n in range(2, m+1)},
          "Missing or out-of-scope graph stratum")
    for allocation in (best, worst[1]):
        _, edges, case = allocation
        check(tutte(edges) == (case["acyclic"], case["cyclic"], case["trees"]),
              "Archived minimum differs under deletion–contraction")
    print(json.dumps(dict(totals), sort_keys=True))
    print(f"Complete strata through {bound} core edges")
    print(f"Minimum scanned ratio: {best[0]}")
    if search is not None:
        print(f"Independent deletion-contraction controls: {controls(search)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ledgers", type=Path, nargs="+")
    parser.add_argument("--search", type=Path)
    parser.add_argument("--geng", type=Path)
    parser.add_argument("--multig", type=Path)
    args = parser.parse_args()
    if bool(args.geng) != bool(args.multig):
        parser.error("provide both geng and multig to verify generator counts")
    programs = tuple(path.resolve(strict=True) for path in (args.geng, args.multig)) if args.geng else None
    verify(args.ledgers, args.search.resolve(strict=True) if args.search else None, programs)
