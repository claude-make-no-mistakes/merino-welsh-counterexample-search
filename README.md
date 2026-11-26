# Multiplicative Merino–Welsh for small multigraphs

**Every loopless bridgeless multigraph with at most 31 edges satisfies the multiplicative Merino–Welsh inequality**

$$T_G(2,0)T_G(0,2)\ \ge\ T_G(1,1)^2.$$

We applied the reduction to 2-inflated matroids in the graphic class. A smallest violation has even size and comes from a core by doubling or subdividing each edge once. At most one edge in each core parallel bundle can be doubled. All eligible decorations on each core are evaluated together by an exact integer tensor transform. [Proof and research note](note.pdf) · [LaTeX](note.tex)

| Core size | Multigraph cores | Covered by the series–parallel theorem | Decorations checked |
| --- | ---: | ---: | ---: |
| [Through 13](through-13.jsonl) | 76,361 | 41,587 | 83,073,664 |
| 14 | 254,671 | 111,233 | 682,969,664 |
| 15 | 1,170,565 | 422,707 | 6,456,573,632 |
| Total | 1,501,597 | 575,527 | 7,222,616,960 |

The computed decorations all passed. Their smallest ratio was $121104/42025$, achieved by subdividing the edges of a spanning path in $K_4$ and doubling the other edges. The ratio describes the computed part. The excluded series–parallel family can attain equality.

## Replay

Use Python 3.10 or later, nauty 2.9.3, and a 64-bit G++ or Clang++ compiler with a 128-bit integer type. Use `search.exe`, `geng.exe`, and `multig.exe` on Windows.

    g++ -std=c++17 -O3 -Wall -Wextra -Werror search.cpp -o search
    python verify.py through-13.jsonl 14-to-15.jsonl --search ./search

The verifier uses deletion–contraction on expanded graphs and independent exact spanning-tree counts. Give it `--geng /path/to/geng --multig /path/to/multig` to reproduce the core count in each stratum. To regenerate the search ledgers from scratch, build `geng` and `multig` from nauty. Use new output paths.

    python enumerate.py replay-13.jsonl --geng /path/to/geng --multig /path/to/multig --search ./search --through 13
    python enumerate.py replay-14-to-15.jsonl --geng /path/to/geng --multig /path/to/multig --search ./search --from 14 --through 15

Nauty generates the simple biconnected supports and then the nonisomorphic positive multiplicity assignments. The separate two-vertex case is handled directly. [Nauty source and documentation](https://pallini.di.uniroma1.it/)

The unrestricted graph conjecture remains open. Knauer, Martínez-Sandoval, and Merino give the reduction to paired expansions and the unbundled rank formulas. Chen, Jin, and Ma give the related minimum-counterexample theorem. The series–parallel theorem handles cores with partial 2-tree supports. Related computations include a large fixed-order, fixed-size multigraph search and Ok's proof for pathwidth at most three. [Reduction and formulas](https://arxiv.org/html/2609.38047v1) · [Minimum-counterexample theorem](https://arxiv.org/html/2610.27706v1#S5) · [Series–parallel theorem](https://arxiv.org/abs/1303.6416) · [Royle's search report](https://symomega.wordpress.com/2013/03/27/more-on-the-merino-welsh-conjecture/) · [Ok's thesis](https://orbit.dtu.dk/en/publications/aspects-of-the-tutte-polynomial/) · [Recent density bounds](https://arxiv.org/html/2611.15008v1)
