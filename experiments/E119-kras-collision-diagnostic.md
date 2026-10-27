# E119: bounded E118 zero-geometry diagnostic

User-returned E118 spread pilot: 10000 conformers, 11.9262 seconds. G12D
2283 and G12V 24547 direct-contact poses; zero poses satisfy the 2.0 A
receptor-heavy-atom distance floor. No soft scores were evaluated. Native self
controls pass, which does not establish library pose-search adequacy.

Question: are failures marginal single-atom contacts or extensive deep overlaps,
and which receptor atoms are involved in the worst pairs?

Replay the exact pilot subset, input seal, pool and finite seeds. Record the
minimum-distance distribution, counts at diagnostic floors 1.0/1.2/1.5/1.8/2.0 A,
number of ligand atoms with nearest-receptor distance below 2.0 A, worst-pair
atom labels, and capped best-boundary examples including coordinates and transforms.
Reconcile both direct and geometry pose counts against E118. Do not change any
acceptance threshold, run a full library scan or export these poses as hits.

These distributions diagnose the distance gate only. Neither failure on all
finite poses nor passing a weakened diagnostic floor establishes binding or
nonbinding. Review marginal examples before selecting refinement or revised
geometry; extensive penetration must not be solved by disabling clashes.
