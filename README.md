# The Solvability Dichotomy for Commutator-Tree Mixing in Finite Groups

Ernest Darell Zermeño — Universidad Panamericana, Campus Guadalajara.

Submitted to the **Journal of the London Mathematical Society**.

This repository contains the manuscript, LaTeX sources and standalone
verification artifact prepared for submission.

The paper identifies when a non-identity iterated commutator forces its
entries to generate a finite group, and derives the resulting mixing law
for random products. For a finite non-solvable group in the configurable
independent-read model, a conjugation-closed generating leaf set has a
uniform strict-contraction level exactly when every proper leaf-generated
subgroup is solvable. The first level is the largest derived length of
those subgroups. The proof combines a fixed-vector criterion with the
normal-generation lemma of Fernández-Alcober and Pintonello.

The companion paper "Uniform Mixing in Barrington Branching Programs" was
published in IACR Communications in Cryptology, Volume 3, Issue 3,
October 2026 (DOI [10.62056/a69qdkmol](https://doi.org/10.62056/a69qdkmol)).

## Paper and reproducibility material

| Location | Contents |
|---|---|
| [paper/main.pdf](paper/main.pdf) | Manuscript (42 pages) |
| [paper/main.tex](paper/main.tex) | LaTeX source (`amsart`), with section files and bibliography |
| [artifact/README.md](artifact/README.md) | Instructions for the 26 verification scripts |
| [artifact/scripts/](artifact/scripts/) | Python scripts: exact certificates, numerical checks, sampled computations |
| [dist/](dist/) | Source ZIP, artifact ZIP, manuscript PDF and SHA-256 checksums |

## Compile the manuscript

With a LaTeX installation providing `pdflatex` and `bibtex`:

```bash
cd paper
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

## Run the verification artifact

Python 3.9+, NumPy 1.24+ and SymPy 1.12+ are required.

```bash
cd artifact
pip install -r requirements.txt
python3 scripts/run_validation.py --output validation-output --jobs 2
```

The runner records every command, exit status, duration and log hash.
Individual scripts can be run from `artifact/scripts/`.

## License

The code is released under the MIT License (see [LICENSE](LICENSE)).
The manuscript text and mathematical content are copyright of the author.
