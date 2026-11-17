from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
import math
import time
from typing import Dict, Mapping, Sequence, Tuple

import numpy as np
import openmm
from openmm import app, unit


ANGSTROM_TO_NM = 0.1
KJ_PER_KCAL = 4.184
DEFAULT_WATER_RESIDUE_NAMES: Tuple[str, ...] = ("HOH", "WAT", "SOL", "TIP3", "TIP3P", "SPC", "SPCE")


@dataclass(frozen=True)
class SCAASWaterMolecule:
    """Index triplet describing one 3-site water molecule in the OpenMM system."""

    oxygen_index: int
    hydrogen1_index: int
    hydrogen2_index: int
    residue_name: str
    residue_id: int


@dataclass(frozen=True)
class SCAASShellDefinition:
    """One polarization shell as used by the Q6 SCAAS algorithm."""

    index: int
    outer_radius_nm: float
    thickness_nm: float
    cstb: float


@dataclass(frozen=True)
class SCAASConfig:
    """Normalized SCAAS parameters stored on EVBSystem after enablement."""

    center_nm: Tuple[float, float, float]
    solvent_radius_nm: float
    temperature_kelvin: float
    shell_radius_nm: float
    shell_force_kj_mol_nm2: float
    excluded_force_kj_mol_nm2: float
    radial_force_kj_mol_nm2: float
    morse_depth_kj_mol: float
    morse_width_per_nm: float
    polarisation_enabled: bool
    charge_correction_enabled: bool
    polarisation_force_kj_mol_rad2: float
    pol_update_steps: int
    water_density_per_a3: float
    water_dipole_ea: float
    dielectric_constant: float
    solvent_residue_names: Tuple[str, ...]
    radial_wall_type: str = "q6_morse"


@dataclass
class SCAASRuntime:
    """Runtime payload needed to attach dynamic SCAAS updates to simulations."""

    config: SCAASConfig
    solvent_molecules: Tuple[SCAASWaterMolecule, ...]
    shells: Tuple[SCAASShellDefinition, ...]
    radial_force: openmm.CustomExternalForce
    shell_fix_force: openmm.CustomExternalForce | None
    polarisation_force: openmm.CustomCompoundBondForce | None
    excluded_atom_indices: Tuple[int, ...]
    shell_atom_indices: Tuple[int, ...]
    frozen_atom_indices: Tuple[int, ...]
    backend_name: str = "cpp"


@dataclass(frozen=True)
class SCAASPolarizationSnapshot:
    """Detailed per-step polarization state for one fixed solvent snapshot."""

    active: Tuple[float, ...]
    targets: Tuple[float, ...]
    corrections: Tuple[float, ...]
    theta_corr: Tuple[float, ...]
    theta_corr_accumulator: Tuple[float, ...]
    theta: Tuple[float, ...]
    radial_nm: Tuple[float, ...]
    shell_members: Tuple[Tuple[int, ...], ...]
    shell_average_theta: Tuple[float | None, ...]
    charge_correction_update: bool


@dataclass(frozen=True)
class Q6RadialBoundaryParameters:
    """Resolved Q6 solvent-boundary parameters in angstrom/kcal units."""

    solvent_radius_angstrom: float
    effective_boundary_angstrom: float
    radial_force_kcal_mol_a2: float
    morse_depth_kcal_mol: float
    morse_width_per_angstrom: float


class SCAASPolarizationReporter:
    """
    Reporter that updates per-water polarization targets every MD step.

    The heavy-force calculation still runs in OpenMM kernels (CustomCompoundBondForce).
    Python updates per-bond shell membership/targets each step, while optional
    charge-correction offsets are refreshed every ``update_interval`` steps.
    """

    def __init__(
        self,
        *,
        force: openmm.CustomCompoundBondForce,
        molecules: Sequence[SCAASWaterMolecule],
        shells: Sequence[SCAASShellDefinition],
        center_nm: Tuple[float, float, float],
        update_interval: int,
        charge_correction: bool,
    ) -> None:
        if update_interval <= 0:
            raise ValueError("'update_interval' must be > 0.")
        self._force = force
        self._molecules = tuple(molecules)
        self._shells = tuple(shells)
        self._center = np.array(center_nm, dtype=float)
        self._update_interval = int(update_interval)
        self._charge_correction = bool(charge_correction)
        self._theta_corr = [0.0 for _ in self._shells]
        self._theta_corr_accumulator = [0.0 for _ in self._shells]
        self._n_updates = 0
        self._last_update_step: int | None = None
        self._n_charge_correction_updates = 0
        self._last_charge_correction_update_step: int | None = None
        self._profile_report_calls = 0
        self._profile_prime_calls = 0
        self._profile_compute_ns = 0
        self._profile_set_bond_ns = 0
        self._profile_update_context_ns = 0
        self._profile_total_report_ns = 0

    @property
    def update_interval(self) -> int:
        return self._update_interval

    @property
    def n_updates(self) -> int:
        return self._n_updates

    @property
    def last_update_step(self) -> int | None:
        return self._last_update_step

    @property
    def theta_corrections(self) -> Tuple[float, ...]:
        return tuple(self._theta_corr)

    @property
    def n_charge_correction_updates(self) -> int:
        return self._n_charge_correction_updates

    @property
    def last_charge_correction_update_step(self) -> int | None:
        return self._last_charge_correction_update_step

    @property
    def profiling_summary(self) -> Mapping[str, object]:
        report_calls = int(self._profile_report_calls)
        mean_divisor = float(report_calls) if report_calls > 0 else 1.0
        return {
            "compute_stage": "python_update",
            "report_calls": report_calls,
            "prime_calls": int(self._profile_prime_calls),
            "compute_ns_total": int(self._profile_compute_ns),
            "set_bond_parameters_ns_total": int(self._profile_set_bond_ns),
            "update_parameters_in_context_ns_total": int(self._profile_update_context_ns),
            "total_report_ns": int(self._profile_total_report_ns),
            "compute_ns_mean": float(self._profile_compute_ns) / mean_divisor if report_calls > 0 else 0.0,
            "set_bond_parameters_ns_mean": float(self._profile_set_bond_ns) / mean_divisor if report_calls > 0 else 0.0,
            "update_parameters_in_context_ns_mean": float(self._profile_update_context_ns) / mean_divisor
            if report_calls > 0
            else 0.0,
            "total_report_ns_mean": float(self._profile_total_report_ns) / mean_divisor if report_calls > 0 else 0.0,
        }

    def describeNextReport(self, simulation: app.Simulation):
        return 1, True, False, False, False, None

    def report(self, simulation: app.Simulation, state: openmm.State) -> None:
        self._update_from_state(simulation, state, update_correction=True)

    def prime(self, simulation: app.Simulation) -> None:
        state = simulation.context.getState(getPositions=True)
        self._update_from_state(simulation, state, update_correction=False)

    def _update_from_state(
        self,
        simulation: app.Simulation,
        state: openmm.State,
        *,
        update_correction: bool,
    ) -> None:
        step = int(simulation.currentStep)
        charge_correction_update = _is_charge_correction_update_step(
            step=step,
            update_correction=bool(update_correction),
            charge_correction=bool(self._charge_correction),
            update_interval=self._update_interval,
            shell_count=len(self._shells),
        )
        total_t0_ns = time.perf_counter_ns()
        positions_nm = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        compute_t0_ns = time.perf_counter_ns()
        active, targets, corrections, theta_corr, theta_corr_accumulator = _compute_polarization_update_python(
            positions_nm=np.asarray(positions_nm, dtype=float),
            molecules=self._molecules,
            shells=self._shells,
            center_nm=tuple(float(x) for x in self._center),
            update_correction=bool(update_correction),
            charge_correction=bool(self._charge_correction),
            step=step,
            update_interval=self._update_interval,
            theta_corr=self._theta_corr,
            theta_corr_accumulator=self._theta_corr_accumulator,
        )
        compute_elapsed_ns = time.perf_counter_ns() - compute_t0_ns
        self._theta_corr = theta_corr
        self._theta_corr_accumulator = theta_corr_accumulator

        set_bond_t0_ns = time.perf_counter_ns()
        for bond_idx, molecule in enumerate(self._molecules):
            self._force.setBondParameters(
                bond_idx,
                (molecule.oxygen_index, molecule.hydrogen1_index, molecule.hydrogen2_index),
                (active[bond_idx], targets[bond_idx], corrections[bond_idx]),
            )
        set_bond_elapsed_ns = time.perf_counter_ns() - set_bond_t0_ns

        update_context_t0_ns = time.perf_counter_ns()
        self._force.updateParametersInContext(simulation.context)
        update_context_elapsed_ns = time.perf_counter_ns() - update_context_t0_ns
        total_elapsed_ns = time.perf_counter_ns() - total_t0_ns

        if update_correction:
            self._n_updates += 1
            self._last_update_step = step
            if charge_correction_update:
                self._n_charge_correction_updates += 1
                self._last_charge_correction_update_step = step
            self._profile_report_calls += 1
            self._profile_compute_ns += int(compute_elapsed_ns)
            self._profile_set_bond_ns += int(set_bond_elapsed_ns)
            self._profile_update_context_ns += int(update_context_elapsed_ns)
            self._profile_total_report_ns += int(total_elapsed_ns)
        else:
            self._profile_prime_calls += 1


def compute_polarization_snapshot_python(
    *,
    positions_nm: np.ndarray,
    molecules: Sequence[SCAASWaterMolecule],
    shells: Sequence[SCAASShellDefinition],
    center_nm: Tuple[float, float, float],
    update_correction: bool,
    charge_correction: bool,
    step: int,
    update_interval: int,
    theta_corr: Sequence[float],
    theta_corr_accumulator: Sequence[float],
) -> SCAASPolarizationSnapshot:
    if update_interval <= 0:
        raise ValueError("'update_interval' must be > 0.")

    theta_corr_local = [float(value) for value in theta_corr]
    theta_corr_accumulator_local = [float(value) for value in theta_corr_accumulator]
    charge_correction_update = _is_charge_correction_update_step(
        step=step,
        update_correction=bool(update_correction),
        charge_correction=bool(charge_correction),
        update_interval=int(update_interval),
        shell_count=len(shells),
    )

    if charge_correction_update:
        for shell_idx, shell in enumerate(shells):
            cstb_clamped = max(-1.0, min(1.0, shell.cstb))
            target_mean = math.acos(cstb_clamped)
            theta_corr_local[shell_idx] += (
                theta_corr_accumulator_local[shell_idx] / float(update_interval) - target_mean
            )
            theta_corr_accumulator_local[shell_idx] = 0.0

    center = np.array(center_nm, dtype=float)
    n_mol = len(molecules)
    theta = [0.0 for _ in range(n_mol)]
    radial = [0.0 for _ in range(n_mol)]

    for idx, molecule in enumerate(molecules):
        o = positions_nm[molecule.oxygen_index]
        h1 = positions_nm[molecule.hydrogen1_index]
        h2 = positions_nm[molecule.hydrogen2_index]
        mu = h1 + h2 - 2.0 * o
        rc = o - center
        mu_norm = float(np.linalg.norm(mu))
        rc_norm = float(np.linalg.norm(rc))
        radial[idx] = rc_norm
        if mu_norm <= 1e-12 or rc_norm <= 1e-12:
            theta[idx] = 0.0
            continue
        scp = float(np.dot(mu, rc) / (mu_norm * rc_norm))
        scp = max(-1.0, min(1.0, scp))
        theta[idx] = math.acos(scp)

    targets = [0.0 for _ in range(n_mol)]
    corrections = [0.0 for _ in range(n_mol)]
    active = [0.0 for _ in range(n_mol)]
    shell_members: list[list[int]] = [[] for _ in shells]
    shell_average_theta: list[float | None] = [None for _ in shells]

    if shells:
        inner_boundary = shells[-1].outer_radius_nm - shells[-1].thickness_nm
        for mol_idx, radius in enumerate(radial):
            if radius <= inner_boundary:
                continue
            shell_idx = 0
            for candidate in range(len(shells) - 1, 0, -1):
                if radius <= shells[candidate].outer_radius_nm:
                    shell_idx = candidate
                    break
            shell_members[shell_idx].append(mol_idx)

    for shell_idx, members in enumerate(shell_members):
        if not members:
            continue
        # Deterministic tie-break rule for equal theta values.
        members.sort(key=lambda idx: (theta[idx], idx))
        n_in_shell = len(members)
        cstb = shells[shell_idx].cstb
        shell_theta_sum = 0.0
        shell_theta_corr = theta_corr_local[shell_idx]
        for rank, mol_idx in enumerate(members, start=1):
            arg = 1.0 + (1.0 - 2.0 * float(rank)) / float(n_in_shell)
            arg = max(-1.0, min(1.0, arg))
            theta_rank = math.acos(arg)
            theta0 = theta_rank - 1.5 * math.sin(theta_rank) * cstb
            theta0 = max(0.0, min(math.pi, theta0))
            targets[mol_idx] = theta0
            corrections[mol_idx] = shell_theta_corr
            active[mol_idx] = 1.0
            shell_theta_sum += theta[mol_idx]

        if update_correction and charge_correction:
            shell_average = shell_theta_sum / float(n_in_shell)
            shell_average_theta[shell_idx] = shell_average
            theta_corr_accumulator_local[shell_idx] += shell_average
        else:
            shell_average_theta[shell_idx] = shell_theta_sum / float(n_in_shell)

    return SCAASPolarizationSnapshot(
        active=tuple(float(value) for value in active),
        targets=tuple(float(value) for value in targets),
        corrections=tuple(float(value) for value in corrections),
        theta_corr=tuple(float(value) for value in theta_corr_local),
        theta_corr_accumulator=tuple(float(value) for value in theta_corr_accumulator_local),
        theta=tuple(float(value) for value in theta),
        radial_nm=tuple(float(value) for value in radial),
        shell_members=tuple(tuple(int(idx) for idx in members) for members in shell_members),
        shell_average_theta=tuple(shell_average_theta),
        charge_correction_update=bool(charge_correction_update),
    )


def _compute_polarization_update_python(
    *,
    positions_nm: np.ndarray,
    molecules: Sequence[SCAASWaterMolecule],
    shells: Sequence[SCAASShellDefinition],
    center_nm: Tuple[float, float, float],
    update_correction: bool,
    charge_correction: bool,
    step: int,
    update_interval: int,
    theta_corr: Sequence[float],
    theta_corr_accumulator: Sequence[float],
) -> Tuple[list[float], list[float], list[float], list[float], list[float]]:
    snapshot = compute_polarization_snapshot_python(
        positions_nm=positions_nm,
        molecules=molecules,
        shells=shells,
        center_nm=center_nm,
        update_correction=update_correction,
        charge_correction=charge_correction,
        step=step,
        update_interval=update_interval,
        theta_corr=theta_corr,
        theta_corr_accumulator=theta_corr_accumulator,
    )
    return (
        list(snapshot.active),
        list(snapshot.targets),
        list(snapshot.corrections),
        list(snapshot.theta_corr),
        list(snapshot.theta_corr_accumulator),
    )


class SCAASCppPolarizationReporter:
    """Reporter using native C++ reference updater for shell/rank/target updates."""

    def __init__(
        self,
        *,
        force: openmm.CustomCompoundBondForce,
        molecules: Sequence[SCAASWaterMolecule],
        shells: Sequence[SCAASShellDefinition],
        center_nm: Tuple[float, float, float],
        update_interval: int,
        charge_correction: bool,
    ) -> None:
        if update_interval <= 0:
            raise ValueError("'update_interval' must be > 0.")
        from . import _scaas_cpp

        self._force = force
        self._molecules = tuple(molecules)
        self._shells = tuple(shells)
        self._center = tuple(float(x) for x in center_nm)
        self._update_interval = int(update_interval)
        self._charge_correction = bool(charge_correction)
        self._theta_corr = np.zeros(len(self._shells), dtype=np.float64)
        self._theta_corr_accumulator = np.zeros(len(self._shells), dtype=np.float64)
        self._n_updates = 0
        self._last_update_step: int | None = None
        self._n_charge_correction_updates = 0
        self._last_charge_correction_update_step: int | None = None
        self._cpp = _scaas_cpp
        bulk_ok, _ = self._cpp.probe_bulk_parameter_updater()
        self._bulk_set_bonds_native = bool(bulk_ok)
        self._profile_report_calls = 0
        self._profile_prime_calls = 0
        self._profile_compute_ns = 0
        self._profile_set_bond_ns = 0
        self._profile_update_context_ns = 0
        self._profile_total_report_ns = 0

        self._oxygen_indices = np.array([mol.oxygen_index for mol in self._molecules], dtype=np.int32)
        self._hydrogen1_indices = np.array([mol.hydrogen1_index for mol in self._molecules], dtype=np.int32)
        self._hydrogen2_indices = np.array([mol.hydrogen2_index for mol in self._molecules], dtype=np.int32)
        self._shell_outer_radius_nm = np.array([shell.outer_radius_nm for shell in self._shells], dtype=np.float64)
        self._shell_thickness_nm = np.array([shell.thickness_nm for shell in self._shells], dtype=np.float64)
        self._shell_cstb = np.array([shell.cstb for shell in self._shells], dtype=np.float64)

    @property
    def update_interval(self) -> int:
        return self._update_interval

    @property
    def n_updates(self) -> int:
        return self._n_updates

    @property
    def last_update_step(self) -> int | None:
        return self._last_update_step

    @property
    def theta_corrections(self) -> Tuple[float, ...]:
        return tuple(float(value) for value in self._theta_corr.tolist())

    @property
    def n_charge_correction_updates(self) -> int:
        return self._n_charge_correction_updates

    @property
    def last_charge_correction_update_step(self) -> int | None:
        return self._last_charge_correction_update_step

    @property
    def profiling_summary(self) -> Mapping[str, object]:
        report_calls = int(self._profile_report_calls)
        mean_divisor = float(report_calls) if report_calls > 0 else 1.0
        return {
            "compute_stage": "cpp_reference_update",
            "report_calls": report_calls,
            "prime_calls": int(self._profile_prime_calls),
            "compute_ns_total": int(self._profile_compute_ns),
            "set_bond_parameters_ns_total": int(self._profile_set_bond_ns),
            "update_parameters_in_context_ns_total": int(self._profile_update_context_ns),
            "total_report_ns": int(self._profile_total_report_ns),
            "compute_ns_mean": float(self._profile_compute_ns) / mean_divisor if report_calls > 0 else 0.0,
            "set_bond_parameters_ns_mean": float(self._profile_set_bond_ns) / mean_divisor if report_calls > 0 else 0.0,
            "update_parameters_in_context_ns_mean": float(self._profile_update_context_ns) / mean_divisor
            if report_calls > 0
            else 0.0,
            "total_report_ns_mean": float(self._profile_total_report_ns) / mean_divisor if report_calls > 0 else 0.0,
        }

    def describeNextReport(self, simulation: app.Simulation):
        return 1, True, False, False, False, None

    def report(self, simulation: app.Simulation, state: openmm.State) -> None:
        self._update_from_state(simulation, state, update_correction=True)

    def prime(self, simulation: app.Simulation) -> None:
        state = simulation.context.getState(getPositions=True)
        self._update_from_state(simulation, state, update_correction=False)

    def _update_from_state(
        self,
        simulation: app.Simulation,
        state: openmm.State,
        *,
        update_correction: bool,
    ) -> None:
        step = int(simulation.currentStep)
        charge_correction_update = _is_charge_correction_update_step(
            step=step,
            update_correction=bool(update_correction),
            charge_correction=bool(self._charge_correction),
            update_interval=self._update_interval,
            shell_count=len(self._shells),
        )
        total_t0_ns = time.perf_counter_ns()
        positions_nm = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        compute_t0_ns = time.perf_counter_ns()
        active, targets, corrections = self._cpp.run_reference_update(
            positions_nm=np.asarray(positions_nm, dtype=np.float64),
            oxygen_indices=self._oxygen_indices,
            hydrogen1_indices=self._hydrogen1_indices,
            hydrogen2_indices=self._hydrogen2_indices,
            shell_outer_radius_nm=self._shell_outer_radius_nm,
            shell_thickness_nm=self._shell_thickness_nm,
            shell_cstb=self._shell_cstb,
            center_nm=self._center,
            update_correction=bool(update_correction),
            charge_correction=self._charge_correction,
            step=step,
            update_interval=self._update_interval,
            theta_corr=self._theta_corr,
            theta_corr_accumulator=self._theta_corr_accumulator,
        )
        compute_elapsed_ns = time.perf_counter_ns() - compute_t0_ns

        set_bond_t0_ns = time.perf_counter_ns()
        if self._bulk_set_bonds_native:
            self._cpp.run_bulk_set_polarization_bonds(
                force=self._force,
                oxygen_indices=self._oxygen_indices,
                hydrogen1_indices=self._hydrogen1_indices,
                hydrogen2_indices=self._hydrogen2_indices,
                active=np.asarray(active, dtype=np.float64),
                targets=np.asarray(targets, dtype=np.float64),
                corrections=np.asarray(corrections, dtype=np.float64),
            )
        else:
            for bond_idx, molecule in enumerate(self._molecules):
                self._force.setBondParameters(
                    bond_idx,
                    (molecule.oxygen_index, molecule.hydrogen1_index, molecule.hydrogen2_index),
                    (float(active[bond_idx]), float(targets[bond_idx]), float(corrections[bond_idx])),
                )
        set_bond_elapsed_ns = time.perf_counter_ns() - set_bond_t0_ns

        update_context_t0_ns = time.perf_counter_ns()
        self._force.updateParametersInContext(simulation.context)
        update_context_elapsed_ns = time.perf_counter_ns() - update_context_t0_ns
        total_elapsed_ns = time.perf_counter_ns() - total_t0_ns

        if update_correction:
            self._n_updates += 1
            self._last_update_step = step
            if charge_correction_update:
                self._n_charge_correction_updates += 1
                self._last_charge_correction_update_step = step
            self._profile_report_calls += 1
            self._profile_compute_ns += int(compute_elapsed_ns)
            self._profile_set_bond_ns += int(set_bond_elapsed_ns)
            self._profile_update_context_ns += int(update_context_elapsed_ns)
            self._profile_total_report_ns += int(total_elapsed_ns)
        else:
            self._profile_prime_calls += 1


@dataclass(frozen=True)
class SCAASBackendResolution:
    """Selected SCAAS backend after applying public selection and availability checks."""

    requested_backend: str
    selected_backend: str
    used_fallback: bool
    fallback_reason: str | None


def is_scaas_cpp_backend_available() -> Tuple[bool, str | None]:
    """Return ``(available, reason)`` for the optional C++ SCAAS backend."""

    try:
        from . import _scaas_cpp
    except ImportError as exc:
        return False, f"module 'evbopenmm._scaas_cpp' is not importable ({exc})"

    ready, detail = _scaas_cpp.probe_reference_updater()
    if not ready:
        return False, f"_scaas_cpp reference updater probe failed ({detail})"
    return True, detail


def resolve_scaas_backend(backend: str) -> SCAASBackendResolution:
    """
    Resolve requested backend string to the supported executable backend.

    Public SCAAS backend selection is now ``cpp`` only.
    """

    requested = str(backend).strip().lower()
    if requested != "cpp":
        raise ValueError(
            "SCAAS backend selection only accepts 'cpp'; "
            f"received '{backend}'. Legacy 'python' and 'auto' backends are no longer public API."
        )

    cpp_available, cpp_reason = is_scaas_cpp_backend_available()
    if not cpp_available:
        reason_text = cpp_reason or "native backend is unavailable"
        raise NotImplementedError(
            "SCAAS backend 'cpp' is not available in this environment. "
            f"Reason: {reason_text}."
        )
    return SCAASBackendResolution(
        requested_backend=requested,
        selected_backend="cpp",
        used_fallback=False,
        fallback_reason=None,
    )


class SCAASBackendController(ABC):
    """Backend-specific runtime attach/diagnostics hooks."""

    @property
    @abstractmethod
    def backend_name(self) -> str:
        """Short backend identifier."""

    @abstractmethod
    def attach(
        self,
        simulation: app.Simulation,
        *,
        runtime: SCAASRuntime,
        evb_scaffold: openmm.System | None,
        initialize: bool,
    ) -> None:
        """Attach backend-specific runtime updates to a simulation."""

    @abstractmethod
    def diagnostics(self, simulation: app.Simulation | None = None) -> Mapping[str, object]:
        """Return backend-specific diagnostics payload."""

    @abstractmethod
    def clear(self) -> None:
        """Clear backend-local runtime state."""


class PythonSCAASBackendController(SCAASBackendController):
    """Current Python reporter-driven SCAAS backend."""

    def __init__(self) -> None:
        self._reporters: Dict[int, SCAASPolarizationReporter] = {}

    @property
    def backend_name(self) -> str:
        return "python"

    def attach(
        self,
        simulation: app.Simulation,
        *,
        runtime: SCAASRuntime,
        evb_scaffold: openmm.System | None,
        initialize: bool,
    ) -> None:
        if simulation.system is not evb_scaffold:
            return
        if _context_parameter_exists(simulation.context, "scaas_enable"):
            simulation.context.setParameter("scaas_enable", 1.0)

        if runtime.polarisation_force is None or not runtime.config.polarisation_enabled:
            return
        key = id(simulation)
        if key in self._reporters:
            if initialize:
                self._reporters[key].prime(simulation)
            return

        reporter = SCAASPolarizationReporter(
            force=runtime.polarisation_force,
            molecules=runtime.solvent_molecules,
            shells=runtime.shells,
            center_nm=runtime.config.center_nm,
            update_interval=runtime.config.pol_update_steps,
            charge_correction=runtime.config.charge_correction_enabled,
        )
        simulation.reporters.append(reporter)
        self._reporters[key] = reporter
        if initialize:
            reporter.prime(simulation)

    def diagnostics(self, simulation: app.Simulation | None = None) -> Mapping[str, object]:
        reporter: SCAASPolarizationReporter | None = None
        if simulation is not None:
            reporter = self._reporters.get(id(simulation))
        elif self._reporters:
            reporter = next(iter(self._reporters.values()))
        if reporter is None:
            return {}
        return {
            "polarisation_updates": reporter.n_updates,
            "polarisation_runtime_updates": reporter.n_updates,
            "last_polarisation_update_step": reporter.last_update_step,
            "charge_correction_updates": reporter.n_charge_correction_updates,
            "last_charge_correction_update_step": reporter.last_charge_correction_update_step,
            "theta_corr": reporter.theta_corrections,
            "profiling": dict(reporter.profiling_summary),
        }

    def clear(self) -> None:
        self._reporters = {}


class CppSCAASBackendController(SCAASBackendController):
    """C++ reference-updater backend controller."""

    def __init__(self) -> None:
        self._reporters: Dict[int, SCAASCppPolarizationReporter] = {}
        self._native_forces: Dict[int, openmm.CustomCompoundBondForce] = {}

    @property
    def backend_name(self) -> str:
        return "cpp"

    def attach(
        self,
        simulation: app.Simulation,
        *,
        runtime: SCAASRuntime,
        evb_scaffold: openmm.System | None,
        initialize: bool,
    ) -> None:
        if simulation.system is not evb_scaffold:
            return
        if _context_parameter_exists(simulation.context, "scaas_enable"):
            simulation.context.setParameter("scaas_enable", 1.0)

        if runtime.polarisation_force is None or not runtime.config.polarisation_enabled:
            return
        key = id(simulation)
        native_force = self._native_forces.get(key)
        if native_force is runtime.polarisation_force:
            return

        native_runtime_attached = False
        if simulation.context.getPlatform().getName().strip().lower() == "cuda":
            from . import _scaas_cpp

            ready, _ = _scaas_cpp.probe_cuda_native_runtime()
            if ready and _scaas_cpp.is_cuda_native_runtime_registered(force=runtime.polarisation_force):
                self._native_forces[key] = runtime.polarisation_force
                self._reporters.pop(key, None)
                native_runtime_attached = True

        if native_runtime_attached:
            return

        self._native_forces.pop(key, None)
        if key in self._reporters:
            if initialize:
                self._reporters[key].prime(simulation)
            return

        reporter = SCAASCppPolarizationReporter(
            force=runtime.polarisation_force,
            molecules=runtime.solvent_molecules,
            shells=runtime.shells,
            center_nm=runtime.config.center_nm,
            update_interval=runtime.config.pol_update_steps,
            charge_correction=runtime.config.charge_correction_enabled,
        )
        simulation.reporters.append(reporter)
        self._reporters[key] = reporter
        if initialize:
            reporter.prime(simulation)

    def diagnostics(self, simulation: app.Simulation | None = None) -> Mapping[str, object]:
        from . import _scaas_cpp

        force: openmm.CustomCompoundBondForce | None = None
        if simulation is not None:
            force = self._native_forces.get(id(simulation))
        elif self._native_forces:
            force = next(iter(self._native_forces.values()))
        if force is not None:
            diagnostics = _scaas_cpp.get_cuda_native_runtime_diagnostics(force=force)
            if diagnostics:
                return diagnostics

        reporter: SCAASCppPolarizationReporter | None = None
        if simulation is not None:
            reporter = self._reporters.get(id(simulation))
        elif self._reporters:
            reporter = next(iter(self._reporters.values()))
        if reporter is None:
            return {}
        return {
            "polarisation_updates": reporter.n_updates,
            "polarisation_runtime_updates": reporter.n_updates,
            "last_polarisation_update_step": reporter.last_update_step,
            "charge_correction_updates": reporter.n_charge_correction_updates,
            "last_charge_correction_update_step": reporter.last_charge_correction_update_step,
            "theta_corr": reporter.theta_corrections,
            "profiling": dict(reporter.profiling_summary),
        }

    def clear(self) -> None:
        self._reporters = {}
        self._native_forces = {}


def build_scaas_backend_controller(backend_name: str) -> SCAASBackendController:
    """Create backend controller by resolved backend name."""

    key = str(backend_name).strip().lower()
    if key == "cpp":
        return CppSCAASBackendController()
    raise ValueError(
        "SCAAS backend controller selection only accepts 'cpp'; "
        f"received '{backend_name}'."
    )


def _context_parameter_exists(context: openmm.Context, parameter_name: str) -> bool:
    try:
        context.getParameter(parameter_name)
    except openmm.OpenMMException:
        return False
    return True


def _is_charge_correction_update_step(
    *,
    step: int,
    update_correction: bool,
    charge_correction: bool,
    update_interval: int,
    shell_count: int,
) -> bool:
    return bool(
        update_correction
        and charge_correction
        and shell_count > 0
        and step > 0
        and update_interval > 0
        and (step % update_interval) == 0
    )


def q6_default_morse_depth(radius_angstrom: float) -> float:
    """Q6 default for solvent Morse depth (kcal/mol)."""

    return 0.26 * math.exp(-0.19 * (float(radius_angstrom) - 15.0)) + 0.74


def q6_default_morse_width(radius_angstrom: float) -> float:
    """Q6 default for solvent Morse width (1/angstrom)."""

    return 0.2 / (1.0 + math.exp(0.4 * (float(radius_angstrom) - 25.0))) + 0.3


def build_q6_radial_boundary_parameters(
    *,
    solvent_radius_angstrom: float,
    temperature_kelvin: float,
    radial_force_kcal_mol_a2: float,
    morse_depth_kcal_mol: float | None = None,
    morse_width_per_angstrom: float | None = None,
) -> Q6RadialBoundaryParameters:
    """
    Resolve the Q6 radial boundary parameters for one solvent radius.

    The effective boundary follows the Q6 manual definition
    ``r0 = rw - sqrt(k_b*T/K)`` in angstrom units.
    """

    solvent_radius_angstrom = float(solvent_radius_angstrom)
    temperature_kelvin = float(temperature_kelvin)
    radial_force_kcal_mol_a2 = float(radial_force_kcal_mol_a2)

    if solvent_radius_angstrom <= 0.0:
        raise ValueError("'solvent_radius_angstrom' must be > 0.")
    if temperature_kelvin < 0.0:
        raise ValueError("'temperature_kelvin' must be >= 0.")
    if radial_force_kcal_mol_a2 <= 0.0:
        raise ValueError("'radial_force_kcal_mol_a2' must be > 0.")

    if morse_depth_kcal_mol is None:
        morse_depth_kcal_mol = q6_default_morse_depth(solvent_radius_angstrom)
    else:
        morse_depth_kcal_mol = float(morse_depth_kcal_mol)
    if morse_width_per_angstrom is None:
        morse_width_per_angstrom = q6_default_morse_width(solvent_radius_angstrom)
    else:
        morse_width_per_angstrom = float(morse_width_per_angstrom)

    kbt_kj_mol = float(
        (unit.MOLAR_GAS_CONSTANT_R * temperature_kelvin * unit.kelvin).value_in_unit(
            unit.kilojoules_per_mole
        )
    )
    shift_angstrom = math.sqrt((kbt_kj_mol / KJ_PER_KCAL) / radial_force_kcal_mol_a2)
    return Q6RadialBoundaryParameters(
        solvent_radius_angstrom=solvent_radius_angstrom,
        effective_boundary_angstrom=solvent_radius_angstrom - shift_angstrom,
        radial_force_kcal_mol_a2=radial_force_kcal_mol_a2,
        morse_depth_kcal_mol=morse_depth_kcal_mol,
        morse_width_per_angstrom=morse_width_per_angstrom,
    )


def q6_radial_boundary_energy_kcal_mol(
    radius_angstrom: float,
    *,
    effective_boundary_angstrom: float,
    radial_force_kcal_mol_a2: float,
    morse_depth_kcal_mol: float,
    morse_width_per_angstrom: float,
) -> float:
    """Evaluate the Q6 radial solvent-boundary potential in kcal/mol."""

    radius_angstrom = float(radius_angstrom)
    effective_boundary_angstrom = float(effective_boundary_angstrom)
    radial_force_kcal_mol_a2 = float(radial_force_kcal_mol_a2)
    morse_depth_kcal_mol = float(morse_depth_kcal_mol)
    morse_width_per_angstrom = float(morse_width_per_angstrom)

    delta = radius_angstrom - effective_boundary_angstrom
    if delta > 0.0:
        return 0.5 * radial_force_kcal_mol_a2 * delta * delta - morse_depth_kcal_mol
    exp_term = math.exp(morse_width_per_angstrom * delta)
    return morse_depth_kcal_mol * (exp_term * exp_term - 2.0 * exp_term)


def build_q6_shells(
    *,
    solvent_radius_angstrom: float,
    total_charge_e: float,
    water_density_per_a3: float,
    water_dipole_ea: float,
    dielectric_constant: float,
    shell_layer_angstrom: float = 3.0001,
    drout_angstrom: float = 0.5,
) -> Tuple[SCAASShellDefinition, ...]:
    """
    Build the outer polarization shells using the same geometry as Q6.

    Returns shells ordered from outermost to innermost.
    """

    if drout_angstrom <= 0.0:
        raise ValueError("'drout_angstrom' must be > 0.")
    if shell_layer_angstrom <= 0.0:
        raise ValueError("'shell_layer_angstrom' must be > 0.")

    drs = float(shell_layer_angstrom) / float(drout_angstrom)
    n_shells = int(-0.5 + math.sqrt(2.0 * drs + 0.25))
    if n_shells <= 0:
        return tuple()

    eps_term = 1.0 - (1.0 / float(dielectric_constant))
    shells: list[SCAASShellDefinition] = []
    rout_a = float(solvent_radius_angstrom)
    for shell_index in range(1, n_shells + 1):
        thickness_a = float(drout_angstrom) * float(shell_index)
        ri_a = rout_a - thickness_a
        rshell_a = (0.5 * (rout_a**3 + ri_a**3)) ** (1.0 / 3.0)
        denominator = float(water_density_per_a3) * float(water_dipole_ea) * 4.0 * math.pi * (rshell_a**2)
        if abs(denominator) <= 1e-12:
            cstb = 0.0
        else:
            cstb = float(total_charge_e) * eps_term / denominator
        shells.append(
            SCAASShellDefinition(
                index=shell_index,
                outer_radius_nm=rout_a * ANGSTROM_TO_NM,
                thickness_nm=thickness_a * ANGSTROM_TO_NM,
                cstb=cstb,
            )
        )
        rout_a = rout_a - thickness_a

    return tuple(shells)


def normalize_scaas_center(
    center: unit.Quantity | Sequence[float] | openmm.Vec3,
) -> Tuple[float, float, float]:
    """Normalize center input to a plain ``(x, y, z)`` tuple in nanometers."""

    if isinstance(center, unit.Quantity):
        values = center.value_in_unit(unit.nanometer)
        x, y, z = float(values[0]), float(values[1]), float(values[2])
        return (x, y, z)

    if isinstance(center, openmm.Vec3):
        return (float(center.x), float(center.y), float(center.z))

    if len(center) != 3:
        raise ValueError("'center' must contain exactly three coordinates.")
    x, y, z = float(center[0]), float(center[1]), float(center[2])
    return (x, y, z)


def normalize_residue_names(names: Sequence[str] | None) -> Tuple[str, ...]:
    if names is None:
        return DEFAULT_WATER_RESIDUE_NAMES
    cleaned: list[str] = []
    for value in names:
        token = str(value).strip().upper()
        if not token:
            continue
        cleaned.append(token)
    if not cleaned:
        raise ValueError("At least one non-empty solvent residue name must be provided.")
    return tuple(dict.fromkeys(cleaned))


def summarize_runtime(
    runtime: SCAASRuntime,
    extra_diagnostics: Mapping[str, object] | None = None,
) -> Mapping[str, object]:
    config = runtime.config
    shell_records = tuple(
        {
            "index": shell.index,
            "outer_radius_nm": shell.outer_radius_nm,
            "thickness_nm": shell.thickness_nm,
            "cstb": shell.cstb,
        }
        for shell in runtime.shells
    )
    payload = {
        "enabled": True,
        "center_nm": config.center_nm,
        "solvent_radius_nm": config.solvent_radius_nm,
        "solvent_molecule_count": len(runtime.solvent_molecules),
        "polarisation_enabled": bool(config.polarisation_enabled),
        "shells": shell_records,
        "excluded_atoms": runtime.excluded_atom_indices,
        "shell_atoms": runtime.shell_atom_indices,
        "backend": runtime.backend_name,
    }
    if extra_diagnostics:
        payload.update(dict(extra_diagnostics))
    return payload
