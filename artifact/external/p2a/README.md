# Historical P2A numerical experiments

The two Python files are copied unchanged from the same author's P2A working tree.
`PROVENANCE.json` records their exact hashes and repository provenance. No separate
LICENSE file was present in the source repository. `s4_family_exact_sup.py`
supplies the exact group routines used by the audit.

From the P1 repository or artifact root, run:

```sh
python3 external/p2a/s4_allpairs_doublecoset_audit.py 11 13 19 29 37
```

The script enumerates double cosets and estimates the largest nontrivial singular
value using 80 power iterations with four starts. These are numerical estimates,
not certified upper bounds. The selected primes reproduce the table in P1.
The source's `generates_full` predicate uses maximality of the exterior S4 for
these primes and a capped closure; it is not a generic generation certificate.
No proof of the all-pairs conjecture is claimed by these experiments.
