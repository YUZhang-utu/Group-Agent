# AIDD Agent end-to-end roadmap

## Mission and operating boundary

The platform is a local-first, resumable design--make--test--learn system. It
starts from immutable molecular and structural evidence, reduces a large
library through validated computational stages, executes approved laboratory
protocols through deterministic adapters, and uses quality-controlled results
to select the next experiment.

The LLM is a planning, explanation, and recommendation layer. It never sends
free-form motion, liquid-handling, temperature, pressure, voltage, or timing
commands directly to a robot or instrument. Hardware executes only a versioned
protocol compiled into a device-specific method after schema validation,
capability checks, safety checks, simulation or dry-run, and the required human
approval.

Experimental raw data and private measurements remain local. Models receive
only explicitly authorized inputs or computed summaries. The separate KRAS
necessity/enhancement work is outside this Project and must not be mixed into
its state.

## System loops

```text
Outer research loop
  evidence synthesis -> hypotheses -> portfolio decision -> new campaign
              ^                                      |
              |                                      v
Inner DMTA loop
  DESIGN -> virtual triage -> MAKE -> TEST -> QC -> LEARN -> next DESIGN
```

Every confirmatory experiment has a locked protocol before execution. Failed
and negative experiments remain first-class results. Each outer-loop review
decides whether to deepen, broaden, pivot, or conclude.

## End-to-end decision graph

```text
Target and campaign definition
  -> receptor/ligand evidence review and immutable locks
  -> library chemistry, microstate, conformer and artifact preparation
  -> broad 2D/3D recall
  -> Gaussian shape/color and reversible anchor reranking
  -> exclusion-volume and limited torsion refinement
  -> docking pose generation and pose QC
  -> independent rescoring and interaction analysis
  -> activity/selectivity/property prediction with uncertainty
  -> molecule-level aggregation, diversity and Pareto selection
  -> human nomination of compounds
  -> inventory/procurement or synthesis planning
  -> executable laboratory protocol and plate map
  -> simulation/dry-run, safety review and human release
  -> robotic execution and instrument acquisition
  -> raw-data registration, assay QC and deterministic analysis
  -> result adjudication and active-learning proposal
  -> human-approved next round
```

No single score is allowed to silently collapse this graph. Every elimination
has a stage, rule, version, reason, and recoverable upstream candidate record.

## Platform planes

### 1. Scientific control plane

- User -> Project -> Task -> Campaign remains the hard ownership boundary.
- A workflow is a versioned DAG of optional Runs, decisions, and Artifacts.
- Human checkpoints lock receptor, pocket, query ligand, docking protocol,
  compound nomination, synthesis order, experimental protocol, and robot run.
- AI recommendations are evidence-bounded and advisory.
- Every Run has an idempotent compute key and explicit lifecycle state.

### 2. Scientific data plane

Maintain immutable, content-addressed artifacts for:

- raw and prepared structures;
- molecule identity, salt/parent relationship, stereochemistry, tautomer and
  protonation microstate;
- conformers, docking poses, trajectories and derived features;
- model inputs, predictions, uncertainty and applicability-domain evidence;
- compounds, batches, samples, containers, wells, transfers and custody;
- protocol source, compiled device methods and execution event logs;
- raw instrument files, parsed measurements, QC decisions and analysis tables.

SQLite stores searchable ownership, identity, state, hashes and lineage. Large
files remain in managed local/project storage or an object store. ELN/LIMS and
inventory systems are adapters; neither becomes an untracked side channel.

### 3. Compute plane

- Workstation: chemistry preparation, visualization, lightweight inference.
- Local SSD: hot immutable library artifacts, indices and reusable caches.
- Slurm/GPU: docking, rescoring, MD and larger model workloads.
- Model registry: versioned activity, selectivity, ADME, safety and pose models.
- Workflow executor: retries, resource limits, manifests and result validation.

### 4. Laboratory execution plane

- Scheduler resolves protocol dependencies, workcell availability, deck
  layout, consumables, calibration status and instrument queues.
- Protocol compiler converts a vendor-neutral assay/synthesis description into
  an allow-listed device method.
- Device adapters control liquid handlers, robotic arms, plate readers, LC-MS,
  incubators and other instruments through supported vendor APIs.
- An event bus records commands, acknowledgements, telemetry, alarms and state
  transitions without treating a command acknowledgment as scientific success.
- A safety supervisor can prevent or abort execution independently of the AI.

### 5. Learning and reporting plane

- Deterministic parsers and QC pipelines transform raw files into measurements.
- Dataset snapshots bind exact molecule/batch/sample/protocol/assay identities.
- Training/validation splits prevent analogue, batch and temporal leakage.
- Models report calibrated uncertainty and applicability domain.
- An acquisition engine proposes the next diverse compound set under potency,
  selectivity, property, cost and information-gain objectives.
- Humans authorize synthesis/procurement and each physical execution round.

## Stage roadmap and gates

### Phase A -- Retrieval foundation (current: E019)

Deliverables:

- accepted real QT9 Reduce angle/projection enrichment;
- Gaussian shape and atom-centered/projected color scoring;
- objective-specific transforms and raw reversible anchor evidence;
- exclusion-volume filtering and terminal-torsion refinement;
- transformed top-hit SDF and complete result manifests;
- Platinum conformer coverage and LIT-PCBA/DUD-E retrieval ablations.

Exit gate:

- multiple WEE1 queries meet preregistered recall and early-enrichment targets;
- no anchor-dependent admission filtering;
- latency, storage, restart and independent-shard append behavior are measured.

### Phase B -- Docking and pose-quality service

Deliverables:

- prepared receptor ensemble with protonation, cofactors, retained waters and
  pocket definitions recorded as immutable artifacts;
- PLANTS adapter first, with a generic interface for later docking engines;
- co-crystal redocking and cross-docking benchmark;
- pose clustering, strain, clash, burial, unsatisfied-polar and interaction-
  fingerprint calculations;
- independent rescoring adapters and molecule-level pose aggregation.

Exit gate:

- redocking/cross-docking acceptance thresholds are locked and met;
- known ligands are enriched over matched controls;
- score direction, failure modes and receptor-ensemble aggregation are explicit;
- docking score is never labeled as experimental affinity.

### Phase C -- Prediction and multi-objective decision layer

Deliverables:

- versioned adapters for activity, kinase selectivity, solubility, permeability,
  clearance, CYP, hERG, Ames, DILI, aggregation and synthesis feasibility;
- standard prediction envelope containing value, units, uncertainty, domain
  status, model/data version and explanation fields;
- temporal/scaffold validation, calibration and model-card registry;
- Pareto ranking, diversity clustering and configurable hard/soft decision rules.

Exit gate:

- every promoted prediction is in-domain or visibly flagged as extrapolation;
- baselines, calibration and prospective evaluation are recorded;
- no unvalidated weighted sum can hide a serious liability;
- users can trace every nominated molecule back to evidence and rejected peers.

### Phase D -- Nomination, inventory and sample identity

Deliverables:

- compound/batch/sample/container/well identity schema;
- barcode support and chain-of-custody event log;
- inventory, concentration, purity, amount, solvent and storage metadata;
- procurement and synthesis-request adapters with explicit human authorization;
- plate-map generator with controls, randomization, replicates and reserved wells.

Exit gate:

- a physical well can be traced to molecule, microstate expectation, batch,
  preparation history and computational nomination;
- unavailable or insufficient material cannot enter an executable protocol;
- plate maps contain positive, negative, blank and technical controls.

### Phase E -- Assay digital twin and manual pilot

Implement one narrow assay before general automation. The recommended first
vertical slice is a low-hazard plate-based binding or biochemical assay already
validated manually in the laboratory.

Deliverables:

- vendor-neutral protocol schema with typed units and allowed parameter ranges;
- equipment/capability registry, labware geometry and deck-layout model;
- protocol linting, volume/mass balance, timing and compatibility checks;
- simulated/dry-run execution and generated operator checklist;
- deterministic raw-data parser and locked QC/analysis protocol;
- manual pilot producing the same artifacts expected from automation.

Exit gate:

- manual reference runs reproduce expected control windows and variability;
- sample/plate/result lineage is complete;
- analysis can distinguish invalid assay runs from inactive compounds;
- no robot is required to prove the information architecture.

### Phase F -- Supervised robotic execution

Automation order:

1. barcode and inventory read-only integration;
2. plate-map and worklist export;
3. liquid-handler dry runs using safe surrogate liquids;
4. plate reader/instrument acquisition and result import;
5. robotic-arm plate transport in a guarded workcell;
6. end-to-end supervised assay execution;
7. limited unattended operation only after reliability evidence.

Deliverables:

- allow-listed device adapters and signed compiled methods;
- preflight checks for calibration, labware, consumables and deck state;
- interlocks, stop/abort handling, recovery SOP and manual takeover;
- command/event/telemetry logs synchronized to a run ID;
- reconciliation of intended versus observed transfers and instrument outputs.

Exit gate:

- repeated supervised runs meet transfer, timing, control-QC and traceability
  targets;
- induced failure tests recover safely without corrupting sample identity;
- a human release token is required for every physical Run;
- network loss, duplicate messages and restart are handled idempotently.

### Phase G -- Closed-loop learning

Deliverables:

- QC-approved assay dataset snapshots and censored/failed-measurement handling;
- active-learning acquisition functions balancing predicted utility, uncertainty,
  novelty, diversity, cost and batch constraints;
- prospective comparison with random, similarity and expert-selection baselines;
- drift, bias, calibration and model-regression monitoring;
- automatic proposal generation, never automatic physical execution.

Exit gate:

- prospective rounds improve a preregistered objective or information-efficiency
  metric over baselines;
- negative and failed results are retained with correct semantics;
- retraining is reproducible from an immutable dataset snapshot;
- the human can reject or edit the proposed next batch without losing lineage.

### Phase H -- Scale, transfer and regulated-readiness options

- multiple assays, workcells, sites and project tenancy;
- role-based access, electronic signatures and audit export where required;
- disaster recovery, backup/restore drills and artifact retention policies;
- equipment qualification, calibration and maintenance evidence;
- change control for models, protocols, adapters and firmware;
- optional integration with organizational ELN/LIMS/QMS systems.

This phase is driven by the intended laboratory and regulatory setting. It must
not be claimed solely from research-prototype tests.

## Required core entities

```text
Molecule -> Microstate -> Conformer -> DockingPose
    |                                      |
    +-> CompoundBatch -> Sample -> Container/Well

ProtocolDefinition -> ProtocolVersion -> CompiledMethod
        -> ExperimentRun -> ExecutionEvent -> RawArtifact
        -> Measurement -> QCDecision -> DatasetSnapshot

ModelVersion -> PredictionRun -> Prediction
DatasetSnapshot -> TrainingRun -> ModelVersion

Campaign -> CandidateDecision -> Nomination -> ExperimentRun
```

Identity must never be inferred from a filename or plate position alone.

## State machines

Computational Run:

```text
draft -> validated -> queued -> running -> succeeded|failed|cancelled
```

Physical Experiment Run:

```text
draft -> protocol_locked -> materials_verified -> simulated
      -> safety_reviewed -> human_released -> executing
      -> completed -> raw_data_registered -> qc_passed|qc_failed
      -> analyzed -> adjudicated
```

Transitions are append-only events. `completed` means equipment execution ended;
it does not mean the assay passed QC or the scientific hypothesis was supported.

## Robot safety and authority boundary

The following are mandatory before connecting physical hardware:

- no direct natural-language-to-motion or natural-language-to-device commands;
- fixed schemas, typed units, parameter ranges and device capability allowlists;
- human approval for materials, protocol version, deck layout and start action;
- vendor safety systems, guards and interlocks remain authoritative;
- independent emergency stop and manual control remain available;
- hazardous chemistry, pressure, heating, biological agents and waste handling
  require institution-approved SOPs and risk assessment outside the Agent;
- simulation and safe-liquid commissioning precede real samples;
- credentials and device-control networks remain outside prompts and source;
- AI cannot suppress alarms, bypass interlocks, widen safety limits or silently
  retry a physical action;
- partial execution is reconciled before any restart to prevent duplicate dosing.

## Metrics dashboard

Computational:

- recall@N, EF1%, BEDROC, diversity, pose RMSD and redocking success;
- model AUROC/PR-AUC where appropriate, RMSE/MAE, calibration, coverage and
  prospective hit rate;
- throughput, latency, storage, cache reuse and failure/restart rate.

Experimental:

- control acceptance, Z-prime or assay-appropriate quality metric, CV, replicate
  agreement, concentration-response fit and invalid-run rate;
- transfer accuracy/precision, sample recovery, contamination/carryover,
  instrument uptime, intervention rate and traceability completeness;
- cycle time, cost per valid measurement and compounds tested per round.

Learning loop:

- valid hits per tested compound and per unit cost/time;
- improvement over random/similarity/expert baselines;
- chemical-space coverage, uncertainty reduction and performance by scaffold;
- reproducibility across batches, plates, days and operators/workcells.

## Near-term implementation order

1. Finish and validate E019; do not dilute the current retrieval milestone.
2. Define E020 docking/redocking protocol and implement the first PLANTS adapter.
3. Define a uniform prediction result envelope and model registry before adding
   many property models.
4. Implement molecule-level evidence aggregation, Pareto selection and a
   nomination checkpoint.
5. Add laboratory identity/provenance entities without controlling hardware.
6. Select one manually validated plate assay and encode its protocol, plate map,
   controls, raw parser and QC rules.
7. Run the same protocol manually and through simulation/dry-run.
8. Integrate one liquid handler and one reader under supervised execution.
9. Close the loop only after prospective assay data and automation QC are sound.

This order prevents unreliable computational rankings or unvalidated assay data
from being amplified by automation.
