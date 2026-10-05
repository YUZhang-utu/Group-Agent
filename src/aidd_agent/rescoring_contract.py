"""Reserved pose-scoring interface; no model loading or execution is registered."""
from dataclasses import dataclass
import math
import re
from typing import Protocol, Sequence


@dataclass(frozen=True)
class PoseInput:
    request_id: str
    molecule_id: str
    conformer_id: str
    receptor_id: str
    receptor_file: str
    receptor_sha256: str
    pose_file: str
    pose_index: int
    pose_file_sha256: str
    microstate_id: str
    source_report_sha256: str


@dataclass(frozen=True)
class ModelSpec:
    name: str
    version: str
    weights_sha256: str
    configuration_sha256: str
    signal: str
    unit: str
    better: str


@dataclass(frozen=True)
class PoseScore:
    request_id: str
    receptor_sha256: str
    pose_file_sha256: str
    pose_index: int
    status: str
    value: float | None
    reason: str = ""


class PoseRescorer(Protocol):
    """Implement once the user's N-E callable and a real I/O example are supplied.

    One adapter invocation produces one signal, with unchanged input coordinates.
    Refinement/cofolding must produce separately tracked structures upstream.
    An adapter must return one result per request, including failures.
    """

    spec: ModelSpec

    def score(self, poses: Sequence[PoseInput]) -> Sequence[PoseScore]: ...


def validate_batch(spec: ModelSpec, poses: Sequence[PoseInput],
                   results: Sequence[PoseScore]) -> dict[str, PoseScore]:
    """Validate attribution and missingness, without certifying chemistry or files.

    The future Project-owned executor must verify file hashes, atom mapping,
    coordinate preservation and model applicability before accepting predictions.
    """
    def digest(value):
        if not isinstance(value, str) or not re.fullmatch(r"[a-f0-9]{64}", value):
            raise ValueError("A lowercase SHA256 digest is required")

    def nonempty(*values):
        if any(not isinstance(v, str) or not v.strip() for v in values):
            raise ValueError("Identifiers and provenance must be nonempty strings")

    nonempty(spec.name, spec.version, spec.unit)
    digest(spec.weights_sha256)
    digest(spec.configuration_sha256)
    if spec.better not in {"higher", "lower"}:
        raise ValueError("Declare score direction")
    if spec.signal not in {"pose_quality", "affinity_prediction", "interaction_score"}:
        raise ValueError("Declare a supported signal family")
    expected = {}
    for pose in poses:
        nonempty(pose.request_id, pose.molecule_id, pose.conformer_id,
                 pose.receptor_id, pose.receptor_file, pose.pose_file, pose.microstate_id)
        for value in (pose.receptor_sha256, pose.pose_file_sha256, pose.source_report_sha256):
            digest(value)
        if type(pose.pose_index) is not int or pose.pose_index < 0:
            raise ValueError("Pose indices must be nonnegative zero-based integers")
        if pose.request_id in expected:
            raise ValueError("Duplicate input request ID")
        expected[pose.request_id] = pose
    checked = {}
    for result in results:
        if result.request_id not in expected or result.request_id in checked:
            raise ValueError("Unknown or duplicate result request ID")
        pose = expected[result.request_id]
        if (type(result.pose_index) is not int or
            (result.receptor_sha256, result.pose_file_sha256, result.pose_index) !=
                (pose.receptor_sha256, pose.pose_file_sha256, pose.pose_index)):
            raise ValueError("Result refers to a different receptor or pose")
        if result.status not in {"ok", "unsupported_input", "invalid_pose", "model_failed", "not_run"}:
            raise ValueError("Unknown scoring status")
        if result.status == "ok":
            if type(result.value) not in {float, int} or not math.isfinite(result.value):
                raise ValueError("Successful scores must be finite numbers")
        elif result.value is not None:
            raise ValueError("Missing or failed scores must be null, never zero-imputed")
        elif not isinstance(result.reason, str) or not result.reason.strip():
            raise ValueError("Unscored inputs need an explicit reason")
        checked[result.request_id] = result
    if checked.keys() != expected.keys():
        raise ValueError("Every input requires a result, including failures")
    return checked
