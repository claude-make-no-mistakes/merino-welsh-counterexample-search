"""Generate all prescribed biconnected multigraph cores using nauty."""

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
from pathlib import Path
import re
import subprocess
import sys


def generate(n, first, last, geng, multig, search):
    if n == 2:
        lines = "".join(f"2 1 0 1 {m}\n" for m in range(first, last + 1)).encode()
        result = subprocess.run([search], input=lines, capture_output=True, check=True)
        expected = last - first + 1
    else:
        with subprocess.Popen([geng, "-q", "-C", str(n), f"{n}:{last}"],
                              stdout=subprocess.PIPE) as simple:
            with subprocess.Popen([multig, "-T", f"-e{first}:{last}"], stdin=simple.stdout,
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE) as weighted:
                simple.stdout.close()
                try:
                    result = subprocess.run([search], stdin=weighted.stdout,
                                            capture_output=True, check=True)
                finally:
                    weighted.stdout.close()
                report = weighted.stderr.read().decode("ascii")
                if weighted.wait() or simple.wait():
                    raise RuntimeError(report)
                found = re.search(r"(\d+) multigraphs written to stdout", report)
                if found is None:
                    raise RuntimeError(report)
                expected = int(found[1])
    rows = [json.loads(line) for line in result.stdout.splitlines()]
    if sum(row["cores"] for row in rows) != expected:
        raise RuntimeError(f"Order {n}: the reader and generator counts disagree")
    print(f"Order {n}: {expected:,} cores", file=sys.stderr, flush=True)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--geng", type=Path, required=True)
    parser.add_argument("--multig", type=Path, required=True)
    parser.add_argument("--search", type=Path, required=True)
    parser.add_argument("--from", type=int, default=2, dest="first")
    parser.add_argument("--through", type=int, required=True, dest="last")
    parser.add_argument("--workers", type=int, default=3)
    args = parser.parse_args()
    if not 2 <= args.first <= args.last <= 20 or args.workers < 1:
        parser.error("expect 2 <= from <= through <= 20 and at least one worker")
    programs = [path.resolve(strict=True) for path in (args.geng, args.multig, args.search)]
    with args.output.open("x", encoding="utf-8", newline="\n") as file:
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            finished = [pool.submit(generate, n, args.first, args.last, *programs)
                        for n in range(2, args.last + 1)]
            rows = []
            for subtask in as_completed(finished):
                rows.extend(subtask.result())
        for row in sorted(rows, key=lambda r: (r["core_edges"], r["vertices"])):
            file.write(json.dumps(row, separators=(",", ":")) + "\n")


if __name__ == "__main__":
    main()
