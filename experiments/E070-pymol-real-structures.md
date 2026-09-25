# E070: deposited chemistry and real-complex interaction replay

Exploratory protocol, 2026-09-25, written before replay.

Hypothesis: relying only on ChemPy bond order 4 loses aromatic assignments when
structures use Kekule bonds. The default also excluded hydrophobic contacts.
Both can leave only yellow polar lines without an execution error.

Compare the existing detector with deposited mmCIF component aromatic flags on
the locally preserved MDM2 structures. Keep coordinates and geometric criteria
unchanged. Report per-complex/per-type counts and source hashes; do not require
every interaction type to occur. Check an independent known positive geometry
and negative controls in unit tests. Check metadata identity mismatches and
removed bonds so source typing cannot override edited chemistry.

Real PyMOL replay should exercise load, selections, detection, distance objects,
colors, export, and repeated invocation. If a compatible runtime is unavailable,
report the limitation and provide a reproducible workstation replay command.
Offline coordinate replay does not validate PyMOL rendering or LLM routing.
