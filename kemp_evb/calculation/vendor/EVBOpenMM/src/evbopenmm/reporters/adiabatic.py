"""AdiabaticReporter — CSV trace of the direct adiabatic-EVB integrator.

Reads the CustomIntegrator's internal globals (E_U, E_S1, E_S2, dV,
c1sq/c2sq, w1/w2, H12, umbrella scalars) at fixed step intervals and
writes one row per interval. The trace is what downstream analysis needs
to reconstruct E_minus(x) plus the mixing-weight time series.

Moved here from ``evbopenmm.sampling.adiabatic`` when the ``reporters``
sub-package was introduced (2026-07-15). ``sampling.adiabatic`` re-exports
``AdiabaticReporter`` from this module for backward compatibility.
"""
from __future__ import annotations

from pathlib import Path
from typing import Tuple

import openmm

_TRACE_FIELDS: Tuple[str, ...] = (
    "step",
    "time_ps",
    "V1_kj_mol",
    "V2_kj_mol",
    "dV_kj_mol",
    "E_minus_kj_mol",
    "c1sq",
    "c2sq",
    "w1",
    "w2",
    "H12_kj_mol",
    "U_umbrella_kj_mol",
    "dV_center_kj_mol",
    "K_umbrella_per_kj_mol",
)


class AdiabaticReporter:
    """OpenMM reporter that writes the adiabatic mixing trace to CSV.

    Reads the integrator's internal globals (``E_U``, ``E_S1``, ``E_S2``,
    ``dV``, ``c1sq``, ``c2sq``, ``w1``, ``w2``, ``H12``) and writes one row
    per ``report_interval_steps``. The first row is preceded by a header.

    The reporter writes the CSV via plain file I/O — no positions or forces
    are pulled off the OpenMM State, so ``describeNextReport`` returns
    ``(steps_until_next, False, False, False, False)``.
    """

    def __init__(
        self,
        file_path,
        *,
        integrator: openmm.CustomIntegrator,
        report_interval_steps: int,
        timestep_ps: float,
        append: bool = False,
    ) -> None:
        if int(report_interval_steps) <= 0:
            raise ValueError("report_interval_steps must be > 0.")
        if float(timestep_ps) <= 0.0:
            raise ValueError("timestep_ps must be > 0.")

        self._integrator = integrator
        self._report_interval_steps = int(report_interval_steps)
        self._timestep_ps = float(timestep_ps)
        path = Path(file_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        self._handle = path.open("a" if append else "w", encoding="utf-8")
        if not append or path.stat().st_size == 0:
            self._handle.write(",".join(_TRACE_FIELDS) + "\n")
            self._handle.flush()

    def describeNextReport(self, simulation):
        current = int(simulation.currentStep)
        steps_until_next = self._report_interval_steps - (
            current % self._report_interval_steps
        )
        return (steps_until_next, False, False, False, False)

    def report(self, simulation, state) -> None:
        e_u = float(self._integrator.getGlobalVariableByName("E_U"))
        e_s1 = float(self._integrator.getGlobalVariableByName("E_S1"))
        e_s2 = float(self._integrator.getGlobalVariableByName("E_S2"))
        dV = float(self._integrator.getGlobalVariableByName("dV"))
        c1sq = float(self._integrator.getGlobalVariableByName("c1sq"))
        c2sq = float(self._integrator.getGlobalVariableByName("c2sq"))
        w1 = float(self._integrator.getGlobalVariableByName("w1"))
        w2 = float(self._integrator.getGlobalVariableByName("w2"))
        h12 = float(self._integrator.getGlobalVariableByName("H12"))
        u_umbrella = float(self._integrator.getGlobalVariableByName("U_umbrella"))
        dv_center = float(self._integrator.getGlobalVariableByName("dV_center"))
        k_umbrella = float(self._integrator.getGlobalVariableByName("K_umbrella"))

        V1 = e_u + e_s1
        V2 = e_u + e_s2
        denom = float(self._integrator.getGlobalVariableByName("denom"))
        alpha = float(self._integrator.getGlobalVariableByName("alpha_kj"))
        e_minus = 0.5 * (V1 + V2 + alpha - denom)

        step = int(simulation.currentStep)
        t_ps = step * self._timestep_ps
        row = (
            step,
            t_ps,
            V1,
            V2,
            dV,
            e_minus,
            c1sq,
            c2sq,
            w1,
            w2,
            h12,
            u_umbrella,
            dv_center,
            k_umbrella,
        )
        self._handle.write(",".join(_format_value(value) for value in row) + "\n")
        self._handle.flush()

    def close(self) -> None:
        if self._handle is not None and not self._handle.closed:
            self._handle.close()

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:
            pass


def _format_value(value) -> str:
    if isinstance(value, int):
        return str(value)
    return f"{float(value):.10g}"


__all__ = ["AdiabaticReporter"]
