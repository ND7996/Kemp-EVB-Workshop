"""EnergyGapReporter — LRA / gap sampling on a running EVB Simulation.

At each ``sample_every_steps`` interval, the reporter:

  1. Calls ``evb_system.set_state_lambdas(sim, probe_hi_lambda)`` and
     reads the potential energy.
  2. Calls ``evb_system.set_state_lambdas(sim, probe_lo_lambda)`` and
     reads it again.
  3. Restores the mapping lambda (``mapping_lambda``) so MD continues
     unaffected.
  4. Optionally repeats steps 1-3 per force group (``per_force_group=True``)
     for a Warshel-style per-group decomposition.
  5. Writes ONE row to a line-buffered CSV.

The CSV is safe to read mid-run (each row is flushed on write; the file
never contains a half-written row). Downstream analysis loads the file
at any point with ``load_gap_samples``.

The three lambda arguments (``mapping_lambda``, ``probe_hi_lambda``,
``probe_lo_lambda``) are forwarded verbatim to ``EVBSystem.set_state_lambdas``.
That method accepts:
    * a scalar float (for the common two-state EVB — sets λ_state1, and
      λ_stateN = 1 − sum),
    * an iterable of N−1 floats in ``EVBSystem.state_names`` order,
    * (test fakes may accept a dict of ``{state_name: weight}``).

ΔU per row = U(probe_hi_lambda) − U(probe_lo_lambda).

For canonical LRA (Warshel):
  * On a *reactant-well* trajectory (mapping λ = ("vb1", 1.0)):
        probe_hi=("vb1",1.0), probe_lo=("vb2",1.0)   →  ⟨ΔU⟩_λ=1
  * On a *transition-window* trajectory (mapping λ = ("vb1",0.5)+("vb2",0.5)):
        same probes → ⟨ΔU⟩_λ=0.5

The reporter never touches positions — the two energy reads happen at
the same instantaneous geometry, only the state weights change between
the two ``sim.context.getState`` calls.

CSV columns (always):

    step, time_ps, U_hi_kj_mol, U_lo_kj_mol, dU_kj_mol

Additional columns if ``per_force_group=True``:

    U_hi_g{k}_kj_mol, U_lo_g{k}_kj_mol, dU_g{k}_kj_mol      for k in 0..31

Only force groups actually assigned to a Force in the System have
meaningful values; unused groups report 0.

Resume: pass ``append=True`` to continue writing to an existing file
(assumes the header matches; caller is responsible for schema stability).
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence, Tuple

import openmm

# ---------------------------------------------------------------------------
# CSV format
# ---------------------------------------------------------------------------

_BASE_FIELDS: Tuple[str, ...] = (
    "step",
    "time_ps",
    "U_hi_kj_mol",
    "U_lo_kj_mol",
    "dU_kj_mol",
)


def _per_group_fields(n_groups: int) -> Tuple[str, ...]:
    fields = []
    for k in range(n_groups):
        fields.extend([
            f"U_hi_g{k}_kj_mol",
            f"U_lo_g{k}_kj_mol",
            f"dU_g{k}_kj_mol",
        ])
    return tuple(fields)


def _format_value(value: Any) -> str:
    if isinstance(value, int):
        return str(value)
    return f"{float(value):.10g}"


# ---------------------------------------------------------------------------
# EVBSystem probing protocol
# ---------------------------------------------------------------------------

def _in_use_force_groups(system: openmm.System) -> Tuple[int, ...]:
    used = set()
    for i in range(system.getNumForces()):
        used.add(int(system.getForce(i).getForceGroup()))
    return tuple(sorted(used))


def _read_group_energies(sim, n_groups: int):
    """Read U per force group by masking groups one at a time. Groups not
    assigned to any Force report 0.

    This does N_groups getState calls per sample; for 5-10 groups the
    per-step overhead is small compared to the MD step at nanosecond scale.
    """
    values = [0.0] * n_groups
    for k in range(n_groups):
        st = sim.context.getState(getEnergy=True, groups={k})
        values[k] = float(st.getPotentialEnergy().value_in_unit(
            openmm.unit.kilojoule_per_mole))
    return values


# ---------------------------------------------------------------------------
# EnergyGapReporter
# ---------------------------------------------------------------------------


class EnergyGapReporter:
    """OpenMM reporter that samples ⟨ΔU⟩ = U(probe_hi) − U(probe_lo).

    Parameters
    ----------
    file_path : str | Path
        CSV file to write. Parent directory is created if missing.
    evb_system : evbopenmm.EVBSystem
        The EVBSystem being propagated. Must expose ``set_state_lambdas(sim, arg)``.
    mapping_lambda : scalar | iterable of float | mapping
        Lambda argument (forwarded to ``evb.set_state_lambdas``) that MD is
        running under. Restored after each sample.
    probe_hi_lambda, probe_lo_lambda : scalar | iterable of float | mapping
        Two lambda combinations to read U at each sample.
    sample_every_steps : int
        MD step cadence between samples.
    timestep_ps : float
        Integrator timestep (used to emit time_ps).
    per_force_group : bool, default False
        If True, also emit per-force-group U_hi / U_lo / dU for each group
        assigned to any Force in the System.
    append : bool, default False
        If True, open the CSV in append mode and skip re-writing the header.
    """

    def __init__(
        self,
        file_path,
        *,
        evb_system,
        mapping_lambda,
        probe_hi_lambda,
        probe_lo_lambda,
        sample_every_steps: int,
        timestep_ps: float,
        per_force_group: bool = False,
        append: bool = False,
    ) -> None:
        if int(sample_every_steps) <= 0:
            raise ValueError("sample_every_steps must be > 0.")
        if float(timestep_ps) <= 0.0:
            raise ValueError("timestep_ps must be > 0.")

        self._evb = evb_system
        self._mapping = mapping_lambda
        self._probe_hi = probe_hi_lambda
        self._probe_lo = probe_lo_lambda
        self._sample_every_steps = int(sample_every_steps)
        self._timestep_ps = float(timestep_ps)
        self._per_force_group = bool(per_force_group)

        # Determine which force groups to emit. We record the fixed set at
        # construction time (from EVBSystem.evb_scaffold) so the CSV schema
        # stays stable across the run.
        self._n_groups = 0
        self._group_ids: Tuple[int, ...] = ()
        if self._per_force_group:
            try:
                scaffold = getattr(evb_system, "evb_scaffold", None)
                if scaffold is not None:
                    self._group_ids = _in_use_force_groups(scaffold)
                    self._n_groups = max(self._group_ids) + 1 if self._group_ids else 0
            except Exception:
                self._n_groups = 0
                self._group_ids = ()

        path = Path(file_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        self._path = path
        self._handle = path.open("a" if append else "w", encoding="utf-8")
        if not append or path.stat().st_size == 0:
            fields = list(_BASE_FIELDS)
            if self._per_force_group:
                fields.extend(_per_group_fields(self._n_groups))
            self._handle.write(",".join(fields) + "\n")
            self._handle.flush()

    # ------------------------------------------------------------------
    # OpenMM Reporter protocol
    # ------------------------------------------------------------------

    def describeNextReport(self, simulation):
        current = int(simulation.currentStep)
        steps_until_next = self._sample_every_steps - (
            current % self._sample_every_steps
        )
        # (steps_until_next, positions, velocities, forces, energies, enforcePeriodicBox)
        # We don't need any State fields — we call getState ourselves after
        # switching λ. Returning False for energies avoids an unnecessary
        # cross-group read.
        return (steps_until_next, False, False, False, False, False)

    def report(self, simulation, state) -> None:
        step = int(simulation.currentStep)
        t_ps = step * self._timestep_ps

        # Probe hi
        self._evb.set_state_lambdas(simulation, self._probe_hi)
        u_hi = float(simulation.context.getState(getEnergy=True)
                     .getPotentialEnergy()
                     .value_in_unit(openmm.unit.kilojoule_per_mole))
        u_hi_groups = _read_group_energies(simulation, self._n_groups) \
            if self._per_force_group else None

        # Probe lo
        self._evb.set_state_lambdas(simulation, self._probe_lo)
        u_lo = float(simulation.context.getState(getEnergy=True)
                     .getPotentialEnergy()
                     .value_in_unit(openmm.unit.kilojoule_per_mole))
        u_lo_groups = _read_group_energies(simulation, self._n_groups) \
            if self._per_force_group else None

        # Restore the mapping λ so MD continues unaffected
        self._evb.set_state_lambdas(simulation, self._mapping)

        du = u_hi - u_lo

        row = [step, t_ps, u_hi, u_lo, du]
        if self._per_force_group:
            for k in range(self._n_groups):
                hi_k = u_hi_groups[k]
                lo_k = u_lo_groups[k]
                row.extend([hi_k, lo_k, hi_k - lo_k])

        self._handle.write(",".join(_format_value(v) for v in row) + "\n")
        self._handle.flush()

    def close(self) -> None:
        if self._handle is not None and not self._handle.closed:
            self._handle.close()

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:
            pass


# ---------------------------------------------------------------------------
# Analysis helpers — kept here so downstream scripts import from the library
# instead of re-implementing block averaging per project.
# ---------------------------------------------------------------------------


def load_gap_samples(path) -> dict:
    """Read a CSV written by EnergyGapReporter into a dict of numpy arrays.

    Returns::

        {
          "step": (N,) int,
          "time_ps": (N,) float,
          "U_hi_kj_mol": (N,) float,
          "U_lo_kj_mol": (N,) float,
          "dU_kj_mol": (N,) float,
          # if per_force_group:
          "U_hi_g0_kj_mol": (N,) float,
          ...
        }
    """
    import numpy as np

    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open(encoding="utf-8") as f:
        header = f.readline().strip().split(",")
    data = np.genfromtxt(
        path, delimiter=",", names=True, skip_header=0, dtype=None, encoding="utf-8",
    )
    # np.genfromtxt with names=True already skips the header row.
    out = {}
    for name in header:
        col = data[name] if data.dtype.names else data
        # ensure numeric dtype; genfromtxt sometimes leaves S dtype for empty col
        if col.dtype.kind in ("S", "U", "O"):
            col = col.astype(float)
        out[name] = col
    # normalise integer column
    if "step" in out:
        out["step"] = out["step"].astype(int)
    return out


def running_mean(values):
    """Cumulative mean of a 1-D array — the running estimate of ⟨X⟩."""
    import numpy as np

    values = np.asarray(values, dtype=float)
    n = len(values)
    if n == 0:
        return values
    idx = np.arange(1, n + 1)
    return np.cumsum(values) / idx


def block_averaged_sem(values, block_sizes: Sequence[int] | None = None):
    """Flyvbjerg-Petersen block-averaging SEM.

    Returns an ``(n_levels, 2)`` ndarray of ``[block_size, SEM]`` rows,
    coarsening the series by successive pair-wise averaging until fewer
    than 4 blocks remain.

    Parameters
    ----------
    values : array-like
        1-D time series (e.g., per-sample ΔU).
    block_sizes : sequence of int, optional
        Ignored — the routine emits the natural coarsening sequence
        (block sizes 1, 2, 4, 8, …). Kept for API symmetry.
    """
    import numpy as np

    x = np.asarray(values, dtype=float).copy()
    if x.size < 4:
        return np.empty((0, 2), dtype=float)
    rows: list[list[float]] = []
    while x.size > 4:
        rows.append([float(x.size), float(x.std(ddof=1) / math.sqrt(x.size))])
        n_pairs = x.size // 2
        x = 0.5 * (x[: 2 * n_pairs : 2] + x[1 : 2 * n_pairs : 2])
    return np.asarray(rows, dtype=float)


__all__ = [
    "EnergyGapReporter",
    "load_gap_samples",
    "running_mean",
    "block_averaged_sem",
]
