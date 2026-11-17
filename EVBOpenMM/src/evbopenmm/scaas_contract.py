from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from typing import Mapping


REFERENCE_PARITY_CONTRACT_NAME = "scaas_reference_parity_v1"
DIAGNOSTICS_PARITY_CONTRACT_NAME = "scaas_diagnostics_parity_v1"
CUDA_LOCKSTEP_TOLERANCE_CONTRACT_NAME = "scaas_cuda_lockstep_v1"
SAMPLED_METRIC_PARITY_CONTRACT_NAME = "scaas_sampled_metric_parity_v1"

# Observed first minimized-start lockstep split under CUDA was 0.01171875 kJ/mol.
# Keep a small margin above that while still treating larger energy-only deltas as regressions.
OBSERVED_CUDA_LOCKSTEP_ENERGY_SPLIT_KJ_MOL = 0.01171875


@dataclass(frozen=True)
class SCAASReferenceParityContract:
    active_abs_tol: float = 1.0e-12
    theta_target_abs_tol: float = 1.0e-12
    theta_correction_abs_tol: float = 1.0e-12
    theta_corr_state_abs_tol: float = 1.0e-12
    theta_accumulator_abs_tol: float = 1.0e-12
    bond_parameter_abs_tol: float = 1.0e-12


@dataclass(frozen=True)
class SCAASDiagnosticsParityContract:
    theta_corr_abs_tol_rad: float = 1.0e-12
    polarisation_updates_must_match: bool = True
    last_update_step_must_match: bool = True


@dataclass(frozen=True)
class SCAASCudaLockstepToleranceContract:
    position_tol_nm: float = 1.0e-7
    velocity_tol_nm_ps: float = 1.0e-7
    force_component_tol_kj_mol_nm: float = 1.0e-7
    theta_tol_rad: float = 1.0e-12
    potential_energy_tol_kj_mol: float = 1.25e-2
    effective_energy_tol_kj_mol: float = 1.25e-2
    reaction_coordinate_tol_kj_mol: float = 1.25e-2
    polarisation_updates_must_match: bool = True
    last_update_step_must_match: bool = True


@dataclass(frozen=True)
class SCAASSampledMetricParityContract:
    energy_abs_tol_kj_mol: float = 1.25e-2
    radius_abs_tol_nm: float = 1.0e-7
    mean_theta_abs_tol_deg: float = 1.0e-7
    shell_count_abs_tol: float = 0.0
    polarisation_updates_must_match: bool = True
    last_update_step_must_match: bool = True


SCAAS_REFERENCE_PARITY_CONTRACT = SCAASReferenceParityContract()
SCAAS_DIAGNOSTICS_PARITY_CONTRACT = SCAASDiagnosticsParityContract()
SCAAS_CUDA_LOCKSTEP_TOLERANCE_CONTRACT = SCAASCudaLockstepToleranceContract()
SCAAS_SAMPLED_METRIC_PARITY_CONTRACT = SCAASSampledMetricParityContract()


def reference_parity_contract_dict() -> dict[str, object]:
    return asdict(SCAAS_REFERENCE_PARITY_CONTRACT)


def diagnostics_parity_contract_dict() -> dict[str, object]:
    return asdict(SCAAS_DIAGNOSTICS_PARITY_CONTRACT)


def cuda_lockstep_tolerance_contract_dict() -> dict[str, object]:
    return asdict(SCAAS_CUDA_LOCKSTEP_TOLERANCE_CONTRACT)


def sampled_metric_parity_contract_dict() -> dict[str, object]:
    return asdict(SCAAS_SAMPLED_METRIC_PARITY_CONTRACT)


def sampled_metric_abs_tolerance(
    metric_name: str,
    contract: SCAASSampledMetricParityContract = SCAAS_SAMPLED_METRIC_PARITY_CONTRACT,
    diagnostics_contract: SCAASDiagnosticsParityContract = SCAAS_DIAGNOSTICS_PARITY_CONTRACT,
) -> float | None:
    name = str(metric_name)
    if name == "polarisation_updates" or name == "last_polarisation_update_step":
        return 0.0
    if name.endswith("_theta_corr_deg"):
        return math.degrees(float(diagnostics_contract.theta_corr_abs_tol_rad))
    if name.endswith("_count") or name == "outer_shell_total_count":
        return float(contract.shell_count_abs_tol)
    if name.endswith("_mean_theta_deg") or name == "theta_all_mean_deg":
        return float(contract.mean_theta_abs_tol_deg)
    if name.endswith("_mean_radius_nm") or name == "radius_all_mean_nm":
        return float(contract.radius_abs_tol_nm)
    if name.endswith("_kj_mol"):
        return float(contract.energy_abs_tol_kj_mol)
    return None


def diagnostics_contract_failures(
    *,
    theta_max_abs_diff: float | None,
    updates_equal: bool,
    last_step_equal: bool,
    contract: SCAASDiagnosticsParityContract = SCAAS_DIAGNOSTICS_PARITY_CONTRACT,
    theta_metric_name: str = "theta_max_abs_diff",
    updates_metric_name: str = "updates_equal",
    last_step_metric_name: str = "last_step_equal",
) -> list[str]:
    failures: list[str] = []

    if theta_max_abs_diff is None or not math.isfinite(float(theta_max_abs_diff)):
        failures.append(theta_metric_name)
    elif float(theta_max_abs_diff) > float(contract.theta_corr_abs_tol_rad):
        failures.append(theta_metric_name)

    if contract.polarisation_updates_must_match and not bool(updates_equal):
        failures.append(updates_metric_name)
    if contract.last_update_step_must_match and not bool(last_step_equal):
        failures.append(last_step_metric_name)
    return failures


def lockstep_contract_failures(
    comparison: Mapping[str, object],
    contract: SCAASCudaLockstepToleranceContract = SCAAS_CUDA_LOCKSTEP_TOLERANCE_CONTRACT,
) -> list[str]:
    reasons: list[str] = []

    def _check(metric: str, tolerance: float) -> None:
        value = comparison.get(metric)
        if value is None:
            return
        if float(value) > float(tolerance):
            reasons.append(metric)

    if not bool(comparison.get("position_shape_equal", False)):
        reasons.append("position_shape_equal")
    if not bool(comparison.get("velocity_shape_equal", False)):
        reasons.append("velocity_shape_equal")
    if not bool(comparison.get("force_shape_equal", False)):
        reasons.append("force_shape_equal")

    _check("position_component_max_abs_diff_nm", float(contract.position_tol_nm))
    _check("velocity_component_max_abs_diff_nm_ps", float(contract.velocity_tol_nm_ps))
    _check("force_component_max_abs_diff_kj_mol_nm", float(contract.force_component_tol_kj_mol_nm))
    _check("theta_corr_max_abs_diff", float(contract.theta_tol_rad))
    _check("potential_energy_abs_diff_kj_mol", float(contract.potential_energy_tol_kj_mol))
    _check("effective_energy_abs_diff_kj_mol", float(contract.effective_energy_tol_kj_mol))
    _check("reaction_coordinate_abs_diff_kj_mol", float(contract.reaction_coordinate_tol_kj_mol))

    if contract.polarisation_updates_must_match and not bool(comparison.get("polarisation_updates_equal", False)):
        reasons.append("polarisation_updates_equal")
    if contract.last_update_step_must_match and not bool(
        comparison.get("last_polarisation_update_step_equal", False)
    ):
        reasons.append("last_polarisation_update_step_equal")
    return reasons


def frame_summary_contract_failures(
    summary: Mapping[str, object],
    contract: SCAASCudaLockstepToleranceContract = SCAAS_CUDA_LOCKSTEP_TOLERANCE_CONTRACT,
) -> list[str]:
    failures: list[str] = []

    def _stat_max(name: str) -> float | None:
        raw = summary.get(name)
        if not isinstance(raw, Mapping):
            return None
        value = raw.get("max")
        if value is None:
            return None
        return float(value)

    def _check(name: str, tolerance: float) -> None:
        value = _stat_max(name)
        if value is None:
            failures.append(name)
            return
        if value > float(tolerance):
            failures.append(name)

    if not bool(summary.get("force_shape_equal_all", False)):
        failures.append("force_shape_equal_all")

    _check("potential_energy_abs_diff_kj_mol", float(contract.potential_energy_tol_kj_mol))
    _check("effective_energy_abs_diff_kj_mol", float(contract.effective_energy_tol_kj_mol))
    _check("reaction_coordinate_abs_diff_kj_mol", float(contract.reaction_coordinate_tol_kj_mol))
    _check("force_component_max_abs_diff_kj_mol_nm", float(contract.force_component_tol_kj_mol_nm))

    failures.extend(
        diagnostics_contract_failures(
            theta_max_abs_diff=_stat_max("theta_corr_max_abs_diff"),
            updates_equal=bool(summary.get("polarisation_updates_equal_all", False)),
            last_step_equal=bool(summary.get("last_polarisation_update_step_equal_all", False)),
            contract=SCAASDiagnosticsParityContract(
                theta_corr_abs_tol_rad=float(contract.theta_tol_rad),
                polarisation_updates_must_match=bool(contract.polarisation_updates_must_match),
                last_update_step_must_match=bool(contract.last_update_step_must_match),
            ),
            theta_metric_name="theta_corr_max_abs_diff",
            updates_metric_name="polarisation_updates_equal_all",
            last_step_metric_name="last_polarisation_update_step_equal_all",
        )
    )
    return failures
