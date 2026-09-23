# E061: avoid unnecessary Gaussian work in contact-first ranking

2026-09-23 protocol before implementation. User reports 2500 of approximately
660000 chunks in ten minutes: naive same-rate projection is 44 hours, not a
measured final runtime. No workstation receipts available yet.

Hypothesis: completion-order scheduling avoids head-of-line stalls; contact-first
ranking needs Gaussian only for poses tied at the maximum contact score within
each conformer. Preserve all seed generation, assignments, pocket checks and
exact contact ties. No candidate, template or conformer budget reductions.

Validate representative equivalence against eager Gaussian, including zero-contact
and tied-contact cases; benchmark a reproducible synthetic pose panel. Record
Gaussian evaluation counts separately. The former all-pose maximum Gaussian
diagnostic becomes explicitly unavailable in lazy mode, not a mislabeled subset
maximum. Production changes require a fresh sealed run; never edit a running
workstation checkout. Workstation speedup requires real receipts and measurements.

Result: 470 passed, 2 skipped. Synthetic scoring panel: 60 shape points drawn
with numpy.default_rng(61), 128 translation seeds, two-feature test geometry;
three warm repetitions. Exact representative equality; Gaussian poses 128 to 1;
median scoring 0.019765s to 0.006963s (2.84x). Excludes seed generation, I/O,
pocket checks and multiprocessing. Artifact: data/e061/synthetic.json.
Workstation receipts and end-to-end speedup remain pending.
