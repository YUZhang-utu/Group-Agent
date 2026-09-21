# E040: installed container AF3 and selectable providers

Implementation acceptance, 2026-09-21. Hypothesis: an explicit Apptainer profile
can preserve the user-provided command semantics while retaining owned task outputs,
compute gates, failed attempts and validated result imports. GPT and DeepSeek must
be independently selectable without putting credentials in repository files.

Validate native compatibility; container argv, binds and XLA; missing host paths;
failed attempts and completed reuse; changed image fingerprint rejection; provider
precedence and credential isolation with mocked HTTP. Run full regression and
English-content/shell checks. No live API, GPU or performance claims from mocks.
Use the existing workstation acceptance procedure for real execution.
