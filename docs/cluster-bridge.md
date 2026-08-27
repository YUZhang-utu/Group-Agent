# Cross-platform Cluster Bridge

The desktop workstation remains the long-running control plane. Cluster access
is performed from an authenticated client session and never requires a daemon,
OpenAI credential, password, MFA token, or private key on a compute node.

## Components

```text
Desktop: aidd-agent control plane and authoritative registry
Client:  OpenSSH on macOS/Linux/Windows or PuTTY on Windows
Cluster: aidd-cluster-runner, registered libraries, environments, and Slurm
```

## Authentication boundary

The bridge generates commands but does not bypass authentication. If a reusable
SSH connection is unavailable, the client must complete the cluster's required
password, MFA, certificate, or Kerberos flow. An already submitted Slurm job
continues after the interactive session expires.

## Cluster Profile

```json
{
  "name": "institutional-hpc",
  "host": "login.hpc.example.edu",
  "port": 22,
  "transport": "auto",
  "ssh_profile": "institutional-hpc",
  "expected_host_key": "SHA256:replace-with-real-fingerprint",
  "project_root": "/project/group/aidd",
  "inbox_root": "/project/group/aidd/inbox",
  "run_root": "/project/group/aidd/runs",
  "library_root": "/project/group/aidd/libraries",
  "software_root": "/project/group/aidd/software",
  "container_root": "/project/group/aidd/containers",
  "scratch_template": "/scratch/{cluster_username}/aidd"
}
```

Passwords and private-key contents are forbidden in this profile.

## Client selection

- macOS: OpenSSH
- Linux: OpenSSH
- Windows: OpenSSH preferred; PuTTY Plink/PSCP supported

OpenSSH uses the registered `ssh_profile` when present. PuTTY requires a saved
session. Connection reuse is opportunistic; failure to reuse a session returns
control to interactive authentication.

## Submission flow

```text
create Run
-> export Slurm Bundle
-> validate Bundle hashes
-> authenticate client
-> upload Bundle
-> execute aidd-cluster-runner submit
-> receive submission_receipt.json
-> import receipt
-> monitor and collect later
```

The Bundle contains structured JSON, a deterministic `submit.slurm`, and a
checksummed manifest. It references registered cluster-side libraries rather
than embedding the full molecular library.

## Security invariants

- Prompt text is never concatenated into a shell command.
- Commands are represented as argument arrays.
- Remote paths must stay inside registered POSIX roots.
- Slurm fields accept restricted character sets and positive resource values.
- A receipt must match the Run, Project, Cluster, and exported Bundle hash.
- Host-key verification must not be disabled.
- Real submission requires an explicit client action after Bundle review.
