from __future__ import annotations

import csv
from collections.abc import Mapping
from dataclasses import asdict, fields, is_dataclass
import json
from pathlib import Path
import shutil
from typing import Any

from openmm import unit

from ..analysis.fitting import (
    fit_evb_parameters_nstate,
    fit_evb_parameters_replicas_nstate,
)
from ..analysis.mapping_fep import analyze_mapping_fep_umbrella_records, write_mapping_fep_outputs
from ..analysis.profile import (
    aggregate_evb_profiles,
    analyze_evb_profile_from_raw,
    analyze_evb_profile_records,
    build_diabatic_profile,
    write_evb_profile_outputs,
    write_evb_replica_outputs,
)
from ..analysis.replicas import analyze_evb_profile_replicas_result
from .layout import ensure_stage_layout, layout_from_analyze_manifest
from .models import (
    AnalyzeWorkflowManifest,
    AnalyzeWorkflowResult,
    AnalyzeWorkflowSpec,
    PrepareWorkflowManifest,
    RunWorkflowManifest,
    StageFileRecord,
)
from .records import read_fep_window_records
from .resume import (
    STAGE_STATUS_COMPLETED,
    STAGE_STATUS_FAILED,
    STAGE_STATUS_RUNNING,
    SubstageStatusRecord,
    stage_resume_state_from_manifest,
    write_stage_resume_state,
)


def analyze_workflow(manifest: AnalyzeWorkflowManifest) -> tuple[dict[str, object], AnalyzeWorkflowManifest]:
    """Execute a Phase B analyze stage from a run-manifest chain."""

    layout = ensure_stage_layout(layout_from_analyze_manifest(manifest))
    current_substage = "load_input_records"
    completed_substages: list[SubstageStatusRecord] = []
    write_stage_resume_state(
        layout,
        stage_resume_state_from_manifest(
            manifest,
            status=STAGE_STATUS_RUNNING,
            active_substage=current_substage,
            substages=(SubstageStatusRecord(current_substage, STAGE_STATUS_RUNNING),),
            detail="Analyze workflow started.",
        ),
    )

    try:
        run_manifest = RunWorkflowManifest.read_json(manifest.spec.run_manifest_path)
        state_names = _resolved_state_names(run_manifest)
        input_records = _resolved_input_record_files(manifest.spec, run_manifest)
        window_records = _load_window_records(input_records)
        completed_substages.append(SubstageStatusRecord(current_substage, STAGE_STATUS_COMPLETED))

        analysis_kind = str(manifest.spec.analysis_kind).strip().lower()
        current_substage = analysis_kind

        if analysis_kind == "mapping_fep_umbrella":
            analysis = _run_mapping_fep_umbrella_analysis(
                manifest=manifest,
                run_manifest=run_manifest,
                state_names=state_names,
                window_records=window_records,
            )
            table_files = _write_mapping_fep_tables(layout, analysis)
            plot_files: tuple[StageFileRecord, ...] = ()
            report_files = (
                _write_analysis_report(
                    layout.report_path("mapping_fep_umbrella_analysis", ".json"),
                    manifest=manifest,
                    run_manifest=run_manifest,
                    input_records=input_records,
                    analysis=analysis,
                ),
            )
        elif analysis_kind == "evb_profile":
            analysis = _run_evb_profile_analysis(
                spec=manifest.spec,
                run_manifest=run_manifest,
                state_names=state_names,
                window_records=window_records,
            )
            mapping_tables = _write_mapping_fep_tables(layout, analysis)
            profile_tables, profile_reports = _write_profile_outputs(layout, analysis)
            table_files = mapping_tables + profile_tables
            plot_files = ()
            report_files = profile_reports + (
                _write_analysis_report(
                    layout.report_path("evb_profile_analysis", ".json"),
                    manifest=manifest,
                    run_manifest=run_manifest,
                    input_records=input_records,
                    analysis=analysis,
                ),
            )
        elif analysis_kind == "evb_profile_replicas":
            replica_window_records = _load_replica_window_records(input_records)
            analysis = _run_evb_profile_replicas_analysis(
                spec=manifest.spec,
                run_manifest=run_manifest,
                state_names=state_names,
                replica_window_records=replica_window_records,
            )
            table_files, replica_reports = _write_replica_outputs(layout, analysis)
            plot_files = ()
            report_files = replica_reports + (
                _write_analysis_report(
                    layout.report_path("evb_profile_replicas_analysis", ".json"),
                    manifest=manifest,
                    run_manifest=run_manifest,
                    input_records=input_records,
                    analysis=analysis,
                ),
            )
        elif analysis_kind == "fit_evb_parameters":
            analysis = _run_fit_evb_parameters_analysis(
                spec=manifest.spec,
                run_manifest=run_manifest,
                state_names=state_names,
                window_records=window_records,
            )
            table_files, fit_reports = _write_fit_outputs(layout, analysis, prefix="fit")
            plot_files = ()
            report_files = fit_reports + (
                _write_analysis_report(
                    layout.report_path("fit_evb_parameters_analysis", ".json"),
                    manifest=manifest,
                    run_manifest=run_manifest,
                    input_records=input_records,
                    analysis=analysis,
                ),
            )
        elif analysis_kind == "fit_evb_parameters_replicas":
            replica_window_records = _load_replica_window_records(input_records)
            analysis = _run_fit_evb_parameters_replicas_analysis(
                spec=manifest.spec,
                run_manifest=run_manifest,
                state_names=state_names,
                replica_window_records=replica_window_records,
            )
            table_files, fit_reports = _write_replica_fit_outputs(layout, analysis)
            plot_files = ()
            report_files = fit_reports + (
                _write_analysis_report(
                    layout.report_path("fit_evb_parameters_replicas_analysis", ".json"),
                    manifest=manifest,
                    run_manifest=run_manifest,
                    input_records=input_records,
                    analysis=analysis,
                ),
            )
        else:
            raise NotImplementedError(
                "Phase B analyze v1 currently supports only "
                "analysis_kind='mapping_fep_umbrella', 'evb_profile', "
                "'evb_profile_replicas', 'fit_evb_parameters', and "
                "'fit_evb_parameters_replicas'."
            )

        analysis_files = tuple(table_files) + tuple(plot_files) + tuple(report_files)
        completed_substages.append(
            SubstageStatusRecord(
                current_substage,
                STAGE_STATUS_COMPLETED,
                generated_files=analysis_files,
            )
        )

        current_substage = "write_analyze_outputs"
        result = AnalyzeWorkflowResult(
            table_files=table_files,
            plot_files=plot_files,
            report_files=report_files,
            manifest_path=layout.manifest_path(),
        )
        output_manifest = AnalyzeWorkflowManifest(
            project=manifest.project,
            spec=manifest.spec,
            result=result,
        )
        output_manifest.write_json(result.manifest_path)
        stage_files = _analyze_stage_output_files(result)
        completed_substages.append(
            SubstageStatusRecord(
                current_substage,
                STAGE_STATUS_COMPLETED,
                generated_files=stage_files,
            )
        )
        _write_completed_stage_resume_state(
            layout,
            output_manifest,
            substages=completed_substages,
            generated_files=stage_files,
        )
        return analysis, output_manifest
    except Exception as exc:  # noqa: BLE001 - record failed stage resume state and re-raise
        _write_failed_stage_resume_state(
            layout,
            manifest,
            current_substage=current_substage,
            completed_substages=completed_substages,
            exc=exc,
        )
        raise


def _run_mapping_fep_umbrella_analysis(
    *,
    manifest: AnalyzeWorkflowManifest,
    run_manifest: RunWorkflowManifest,
    state_names: tuple[str, ...],
    window_records,
) -> dict[str, object]:
    return analyze_mapping_fep_umbrella_records(
        window_records,
        state_names=state_names,
        kT_kj_mol=_thermal_energy_kj_mol(manifest.spec, run_manifest),
        skip_samples=int(manifest.spec.skip_samples),
        n_gap_bins=int(manifest.spec.n_gap_bins or 40),
        min_points_per_bin=int(manifest.spec.min_points_per_bin or 20),
        input_energy_unit="kJ/mol",
    )


def _run_evb_profile_analysis(
    *,
    spec: AnalyzeWorkflowSpec,
    run_manifest: RunWorkflowManifest,
    state_names: tuple[str, ...],
    window_records,
) -> dict[str, object]:
    base_analysis = analyze_mapping_fep_umbrella_records(
        window_records,
        state_names=state_names,
        kT_kj_mol=_thermal_energy_kj_mol(spec, run_manifest),
        skip_samples=int(spec.skip_samples),
        n_gap_bins=int(spec.n_gap_bins or 40),
        min_points_per_bin=int(spec.min_points_per_bin or 20),
        input_energy_unit="kJ/mol",
    )
    if len(state_names) != 2:
        raise ValueError(
            "analysis_kind='evb_profile' requires exactly two diabatic states."
        )
    alpha = float(spec.alpha_kcal_mol or 0.0)
    hij = abs(float(spec.hij_kcal_mol or 0.0))
    return analyze_evb_profile_records(
        window_records,
        state_names=state_names,
        base_analysis=base_analysis,
        diabatic_offsets_kcal_mol=(0.0, alpha),
        coupling_matrix_kcal_mol=((0.0, hij), (hij, 0.0)),
        reaction_coordinate_states=(0, 1),
        skip_samples=int(spec.skip_samples),
        n_gap_bins=int(spec.n_gap_bins or 40),
        min_points_per_bin=int(spec.min_points_per_bin or 20),
        input_energy_unit="kJ/mol",
    )


def _run_evb_profile_replicas_analysis(
    *,
    spec: AnalyzeWorkflowSpec,
    run_manifest: RunWorkflowManifest,
    state_names: tuple[str, ...],
    replica_window_records: tuple[tuple[object, ...], ...],
) -> dict[str, object]:
    return analyze_evb_profile_replicas_result(
        replica_window_records,
        analyze_profile=lambda records: _run_evb_profile_analysis(
            spec=spec,
            run_manifest=run_manifest,
            state_names=state_names,
            window_records=records,
        ),
        aggregate_profiles=aggregate_evb_profiles,
        minimum_replicas_per_bin=max(int(spec.minimum_replicas_per_bin), 1),
        output_directory=None,
    )


def _state_index_map(state_names: tuple[str, ...]) -> dict[str, int]:
    return {name: idx for idx, name in enumerate(state_names)}


def _assemble_diabatic_parameters(
    state_names: tuple[str, ...],
    couplings: Mapping[tuple[str, str], float],
    offsets: Mapping[str, float],
) -> tuple[tuple[float, ...], tuple[tuple[float, ...], ...]]:
    index = _state_index_map(state_names)
    n_states = len(state_names)
    offset_vec = [0.0] * n_states
    for name, value in offsets.items():
        offset_vec[index[name]] = float(value)
    coupling_mat = [[0.0] * n_states for _ in range(n_states)]
    for (i, j), value in couplings.items():
        a = index[i]
        b = index[j]
        coupling_mat[a][b] = float(value)
        coupling_mat[b][a] = float(value)
    return tuple(offset_vec), tuple(tuple(row) for row in coupling_mat)


def _resolve_fit_workflow_defaults(
    spec: AnalyzeWorkflowSpec,
    state_names: tuple[str, ...],
) -> tuple[tuple[str, ...], tuple[int, int]]:
    offset_pattern = (
        tuple(spec.diabatic_offset_pattern)
        if spec.diabatic_offset_pattern is not None
        else tuple(state_names[1:])
    )
    reaction_coordinate_states = (
        tuple(spec.reaction_coordinate_states)
        if spec.reaction_coordinate_states is not None
        else (0, len(state_names) - 1)
    )
    return offset_pattern, reaction_coordinate_states


def _run_fit_evb_parameters_analysis(
    *,
    spec: AnalyzeWorkflowSpec,
    run_manifest: RunWorkflowManifest,
    state_names: tuple[str, ...],
    window_records,
) -> dict[str, object]:
    kT_kj = _thermal_energy_kj_mol(spec, run_manifest)
    n_gap_bins = int(spec.n_gap_bins or 40)
    min_points_per_bin = int(spec.min_points_per_bin or 20)
    skip_samples = int(spec.skip_samples)
    base_analysis = analyze_mapping_fep_umbrella_records(
        window_records,
        state_names=state_names,
        kT_kj_mol=kT_kj,
        skip_samples=skip_samples,
        n_gap_bins=n_gap_bins,
        min_points_per_bin=min_points_per_bin,
        input_energy_unit="kJ/mol",
    )
    offset_pattern, rcs = _resolve_fit_workflow_defaults(spec, state_names)

    # The raw diabatic profile is alpha/H12-INDEPENDENT, so bin it ONCE here and
    # re-run only the cheap analyser (analyze_evb_profile_from_raw) per optimiser trial.
    raw_profile = build_diabatic_profile(
        window_records,
        state_names=state_names,
        base_analysis=base_analysis,
        reaction_coordinate_states=rcs,
        skip_samples=skip_samples,
        n_gap_bins=n_gap_bins,
        min_points_per_bin=min_points_per_bin,
        input_energy_unit="kJ/mol",
    )

    def _analyze_candidate(
        couplings: Mapping[tuple[str, str], float],
        offsets: Mapping[str, float],
    ) -> Mapping[str, object]:
        offset_vec, coupling_mat = _assemble_diabatic_parameters(state_names, couplings, offsets)
        return analyze_evb_profile_from_raw(
            raw_profile,
            window_records=window_records,
            base_analysis=base_analysis,
            state_names=state_names,
            skip_samples=skip_samples,
            diabatic_offsets_kcal_mol=offset_vec,
            coupling_matrix_kcal_mol=coupling_mat,
        )

    return fit_evb_parameters_nstate(
        state_names=state_names,
        coupling_pattern=spec.coupling_pattern,
        diabatic_offset_pattern=offset_pattern,
        analyze_candidate=_analyze_candidate,
        target_dg_dagger_kcal_mol=spec.target_dg_dagger_kcal_mol,
        target_dg0_kcal_mol=spec.target_dg0_kcal_mol,
        weight_dg_dagger=float(spec.weight_dg_dagger),
        weight_dg0=float(spec.weight_dg0),
        coupling_bounds_kcal_mol=spec.coupling_bounds_kcal_mol,
        offset_bounds_kcal_mol=spec.offset_bounds_kcal_mol,
        scipy_method=str(spec.scipy_method),
        scipy_tolerance=spec.scipy_tolerance,
    )


def _run_fit_evb_parameters_replicas_analysis(
    *,
    spec: AnalyzeWorkflowSpec,
    run_manifest: RunWorkflowManifest,
    state_names: tuple[str, ...],
    replica_window_records: tuple[tuple[object, ...], ...],
) -> dict[str, object]:
    kT_kj = _thermal_energy_kj_mol(spec, run_manifest)
    n_gap_bins = int(spec.n_gap_bins or 40)
    min_points_per_bin = int(spec.min_points_per_bin or 20)
    skip_samples = int(spec.skip_samples)
    minimum_replicas_per_bin = max(int(spec.minimum_replicas_per_bin), 1)
    offset_pattern, rcs = _resolve_fit_workflow_defaults(spec, state_names)

    # Build each replica's alpha/H12-INDEPENDENT base + raw profile ONCE; per optimiser
    # trial re-run only the cheap analyser on the cached raw (keyed by records identity).
    replica_cache: dict[int, tuple[Mapping[str, object], Mapping[str, object]]] = {}
    for records in replica_window_records:
        base = analyze_mapping_fep_umbrella_records(
            records,
            state_names=state_names,
            kT_kj_mol=kT_kj,
            skip_samples=skip_samples,
            n_gap_bins=n_gap_bins,
            min_points_per_bin=min_points_per_bin,
            input_energy_unit="kJ/mol",
        )
        raw = build_diabatic_profile(
            records,
            state_names=state_names,
            base_analysis=base,
            reaction_coordinate_states=rcs,
            skip_samples=skip_samples,
            n_gap_bins=n_gap_bins,
            min_points_per_bin=min_points_per_bin,
            input_energy_unit="kJ/mol",
        )
        replica_cache[id(records)] = (base, raw)

    def _analyze_candidate(
        couplings: Mapping[tuple[str, str], float],
        offsets: Mapping[str, float],
    ) -> Mapping[str, object]:
        offset_vec, coupling_mat = _assemble_diabatic_parameters(state_names, couplings, offsets)

        def _analyze_one(records) -> Mapping[str, object]:
            base, raw = replica_cache[id(records)]
            return analyze_evb_profile_from_raw(
                raw,
                window_records=records,
                base_analysis=base,
                state_names=state_names,
                diabatic_offsets_kcal_mol=offset_vec,
                coupling_matrix_kcal_mol=coupling_mat,
                skip_samples=skip_samples,
            )

        return analyze_evb_profile_replicas_result(
            replica_window_records,
            analyze_profile=_analyze_one,
            aggregate_profiles=aggregate_evb_profiles,
            minimum_replicas_per_bin=minimum_replicas_per_bin,
            output_directory=None,
        )

    return fit_evb_parameters_replicas_nstate(
        state_names=state_names,
        replica_count=len(replica_window_records),
        coupling_pattern=spec.coupling_pattern,
        diabatic_offset_pattern=offset_pattern,
        analyze_candidate=_analyze_candidate,
        target_dg_dagger_kcal_mol=spec.target_dg_dagger_kcal_mol,
        target_dg0_kcal_mol=spec.target_dg0_kcal_mol,
        weight_dg_dagger=float(spec.weight_dg_dagger),
        weight_dg0=float(spec.weight_dg0),
        coupling_bounds_kcal_mol=spec.coupling_bounds_kcal_mol,
        offset_bounds_kcal_mol=spec.offset_bounds_kcal_mol,
        scipy_method=str(spec.scipy_method),
        scipy_tolerance=spec.scipy_tolerance,
    )


def _write_mapping_fep_tables(layout, analysis: Mapping[str, object]) -> tuple[StageFileRecord, ...]:
    write_mapping_fep_outputs(analysis, output_directory=layout.table_artifacts_dir)
    return (
        StageFileRecord("fep_windows_csv", layout.table_artifacts_dir / "fep_windows.csv"),
        StageFileRecord("fep_edges_csv", layout.table_artifacts_dir / "fep_edges.csv"),
        StageFileRecord("umbrella_profile_csv", layout.table_artifacts_dir / "umbrella_profile.csv"),
    )


def _write_profile_outputs(layout, analysis: Mapping[str, object]) -> tuple[tuple[StageFileRecord, ...], tuple[StageFileRecord, ...]]:
    staging_dir = layout.artifacts_dir / "_profile_staging"
    if staging_dir.exists():
        shutil.rmtree(staging_dir)
    staging_dir.mkdir(parents=True, exist_ok=True)
    try:
        write_evb_profile_outputs(analysis, output_directory=staging_dir)
        profile_table_source = staging_dir / "evb_profile.csv"
        thermo_report_source = staging_dir / "evb_thermochemistry.txt"

        profile_table_path = layout.table_artifacts_dir / profile_table_source.name
        thermo_report_path = layout.report_artifacts_dir / thermo_report_source.name
        profile_table_path.parent.mkdir(parents=True, exist_ok=True)
        thermo_report_path.parent.mkdir(parents=True, exist_ok=True)

        shutil.move(str(profile_table_source), profile_table_path)
        shutil.move(str(thermo_report_source), thermo_report_path)
    finally:
        if staging_dir.exists():
            shutil.rmtree(staging_dir)

    return (
        (StageFileRecord("evb_profile_csv", profile_table_path),),
        (StageFileRecord("evb_thermochemistry_txt", thermo_report_path),),
    )


def _write_replica_outputs(layout, analysis: Mapping[str, object]) -> tuple[tuple[StageFileRecord, ...], tuple[StageFileRecord, ...]]:
    replica_dir = layout.artifacts_dir / "replicas"
    if replica_dir.exists():
        shutil.rmtree(replica_dir)
    replica_dir.mkdir(parents=True, exist_ok=True)
    for idx, replica_analysis in enumerate(analysis["replicas"]):
        replica_out = replica_dir / f"rep_{idx:03d}"
        write_mapping_fep_outputs(replica_analysis, output_directory=replica_out)
        write_evb_profile_outputs(replica_analysis, output_directory=replica_out)

    staging_dir = layout.artifacts_dir / "_replica_aggregate_staging"
    if staging_dir.exists():
        shutil.rmtree(staging_dir)
    staging_dir.mkdir(parents=True, exist_ok=True)
    try:
        write_evb_replica_outputs(analysis["aggregate"], output_directory=staging_dir)
        ensemble_source = staging_dir / "evb_replica_profile_ensemble.csv"
        thermo_source = staging_dir / "evb_replica_thermochemistry.txt"
        ensemble_path = layout.table_artifacts_dir / ensemble_source.name
        thermo_path = layout.report_artifacts_dir / thermo_source.name
        ensemble_path.parent.mkdir(parents=True, exist_ok=True)
        thermo_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(ensemble_source), ensemble_path)
        shutil.move(str(thermo_source), thermo_path)
    finally:
        if staging_dir.exists():
            shutil.rmtree(staging_dir)

    return (
        (StageFileRecord("evb_replica_profile_ensemble_csv", ensemble_path),),
        (StageFileRecord("evb_replica_thermochemistry_txt", thermo_path),),
    )


def _write_fit_outputs(
    layout,
    fit_result: Mapping[str, object],
    *,
    prefix: str,
) -> tuple[tuple[StageFileRecord, ...], tuple[StageFileRecord, ...]]:
    history_path = layout.table_artifacts_dir / f"{prefix}_history.csv"
    _write_fit_history_csv(history_path, fit_result)
    best_tables, best_reports = _write_best_profile_analysis_outputs(
        layout,
        fit_result["best_analysis"],
        prefix=f"{prefix}_best",
    )
    return (
        (StageFileRecord(f"{prefix}_history_csv", history_path),) + best_tables,
        best_reports,
    )


def _write_replica_fit_outputs(
    layout,
    fit_result: Mapping[str, object],
) -> tuple[tuple[StageFileRecord, ...], tuple[StageFileRecord, ...]]:
    history_path = layout.table_artifacts_dir / "replica_fit_history.csv"
    _write_fit_history_csv(history_path, fit_result)
    best_tables, best_reports = _write_best_replica_analysis_outputs(
        layout,
        fit_result["best_replica_result"],
        prefix="replica_fit_best",
    )
    return (
        (StageFileRecord("replica_fit_history_csv", history_path),) + best_tables,
        best_reports,
    )


def _write_best_profile_analysis_outputs(
    layout,
    analysis: Mapping[str, object],
    *,
    prefix: str,
) -> tuple[tuple[StageFileRecord, ...], tuple[StageFileRecord, ...]]:
    staging_dir = layout.artifacts_dir / f"_{prefix}_profile_staging"
    if staging_dir.exists():
        shutil.rmtree(staging_dir)
    staging_dir.mkdir(parents=True, exist_ok=True)
    try:
        write_mapping_fep_outputs(analysis, output_directory=staging_dir)
        write_evb_profile_outputs(analysis, output_directory=staging_dir)
        table_map = {
            f"{prefix}_fep_windows_csv": (staging_dir / "fep_windows.csv", layout.table_artifacts_dir / f"{prefix}_fep_windows.csv"),
            f"{prefix}_fep_edges_csv": (staging_dir / "fep_edges.csv", layout.table_artifacts_dir / f"{prefix}_fep_edges.csv"),
            f"{prefix}_umbrella_profile_csv": (staging_dir / "umbrella_profile.csv", layout.table_artifacts_dir / f"{prefix}_umbrella_profile.csv"),
            f"{prefix}_evb_profile_csv": (staging_dir / "evb_profile.csv", layout.table_artifacts_dir / f"{prefix}_evb_profile.csv"),
        }
        report_map = {
            f"{prefix}_evb_thermochemistry_txt": (
                staging_dir / "evb_thermochemistry.txt",
                layout.report_artifacts_dir / f"{prefix}_evb_thermochemistry.txt",
            ),
        }
        table_files = _move_stage_outputs(table_map)
        report_files = _move_stage_outputs(report_map)
    finally:
        if staging_dir.exists():
            shutil.rmtree(staging_dir)

    return table_files, report_files


def _write_best_replica_analysis_outputs(
    layout,
    replica_result: Mapping[str, object],
    *,
    prefix: str,
) -> tuple[tuple[StageFileRecord, ...], tuple[StageFileRecord, ...]]:
    staging_dir = layout.artifacts_dir / f"_{prefix}_replica_staging"
    if staging_dir.exists():
        shutil.rmtree(staging_dir)
    staging_dir.mkdir(parents=True, exist_ok=True)
    try:
        write_evb_replica_outputs(replica_result["aggregate"], output_directory=staging_dir)
        table_map = {
            f"{prefix}_profile_ensemble_csv": (
                staging_dir / "evb_replica_profile_ensemble.csv",
                layout.table_artifacts_dir / f"{prefix}_profile_ensemble.csv",
            ),
        }
        report_map = {
            f"{prefix}_thermochemistry_txt": (
                staging_dir / "evb_replica_thermochemistry.txt",
                layout.report_artifacts_dir / f"{prefix}_thermochemistry.txt",
            ),
        }
        table_files = _move_stage_outputs(table_map)
        report_files = _move_stage_outputs(report_map)
    finally:
        if staging_dir.exists():
            shutil.rmtree(staging_dir)

    return table_files, report_files


def _move_stage_outputs(output_map: Mapping[str, tuple[Path, Path]]) -> tuple[StageFileRecord, ...]:
    records: list[StageFileRecord] = []
    for label, (source, destination) in output_map.items():
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(source), destination)
        records.append(StageFileRecord(label, destination))
    return tuple(records)


def _coupling_column_name(pair: tuple[str, str]) -> str:
    a, b = pair
    return f"Hij_{a}__{b}_kcal_mol"


def _offset_column_name(state_name: str) -> str:
    return f"offset_{state_name}_kcal_mol"


def _write_fit_history_csv(
    path: Path,
    fit_result: Mapping[str, object],
) -> None:
    coupling_pattern = tuple(tuple(pair) for pair in fit_result["coupling_pattern"])
    offset_pattern = tuple(fit_result["diabatic_offset_pattern"])
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        header: list[str] = ["iteration", "objective_kcal2_mol2", "dg_dagger_kcal_mol", "dg0_kcal_mol"]
        for pair in coupling_pattern:
            header.append(_coupling_column_name(pair))
        for state_name in offset_pattern:
            header.append(_offset_column_name(state_name))
        writer.writerow(header)

        for idx, entry in enumerate(fit_result["history"], start=1):
            row: list[object] = [
                idx,
                _fmt_optional(entry.get("objective_kcal2_mol2")),
                _fmt_optional(entry.get("dg_dagger_kcal_mol")),
                _fmt_optional(entry.get("dg0_kcal_mol")),
            ]
            couplings = entry.get("couplings_kcal_mol", {})
            offsets = entry.get("offsets_kcal_mol", {})
            for pair in coupling_pattern:
                value = couplings.get(pair)
                if value is None:
                    value = couplings.get((pair[1], pair[0]))
                row.append(_fmt_optional(value))
            for state_name in offset_pattern:
                row.append(_fmt_optional(offsets.get(state_name)))
            writer.writerow(row)


def _write_analysis_report(
    path: Path,
    *,
    manifest: AnalyzeWorkflowManifest,
    run_manifest: RunWorkflowManifest,
    input_records: tuple[StageFileRecord, ...],
    analysis: Mapping[str, object],
) -> StageFileRecord:
    report_payload = {
        "analysis_kind": manifest.spec.analysis_kind,
        "run_manifest_path": str(manifest.spec.run_manifest_path),
        "input_record_files": [record.to_dict() for record in input_records],
        "state_names": list(_resolved_state_names(run_manifest)),
        "analysis": analysis,
    }
    _write_json(path, report_payload)
    return StageFileRecord(f"{path.stem}_json", path)


def _write_completed_stage_resume_state(
    layout,
    manifest: AnalyzeWorkflowManifest,
    *,
    substages: list[SubstageStatusRecord],
    generated_files: tuple[StageFileRecord, ...],
) -> Path:
    return write_stage_resume_state(
        layout,
        stage_resume_state_from_manifest(
            manifest,
            status=STAGE_STATUS_COMPLETED,
            substages=tuple(substages),
            generated_files=generated_files,
            detail="Analyze workflow completed successfully.",
        ),
    )


def _write_failed_stage_resume_state(
    layout,
    manifest: AnalyzeWorkflowManifest,
    *,
    current_substage: str,
    completed_substages: list[SubstageStatusRecord],
    exc: Exception,
) -> Path:
    failed_substages = tuple(completed_substages) + (SubstageStatusRecord(current_substage, STAGE_STATUS_FAILED),)
    return write_stage_resume_state(
        layout,
        stage_resume_state_from_manifest(
            manifest,
            status=STAGE_STATUS_FAILED,
            active_substage=current_substage,
            substages=failed_substages,
            detail="Analyze workflow failed before completion.",
            error_message=f"{type(exc).__name__}: {exc}",
        ),
    )


def _analyze_stage_output_files(result: AnalyzeWorkflowResult) -> tuple[StageFileRecord, ...]:
    records: list[StageFileRecord] = []
    if result.manifest_path is not None:
        records.append(StageFileRecord("analyze_manifest_json", result.manifest_path))
    records.extend(result.table_files)
    records.extend(result.plot_files)
    records.extend(result.report_files)
    return tuple(records)


def _resolved_state_names(run_manifest: RunWorkflowManifest) -> tuple[str, ...]:
    if run_manifest.spec.states:
        return tuple(run_manifest.spec.states)

    prepare_manifest = PrepareWorkflowManifest.read_json(run_manifest.spec.prepared_manifest_path)
    if prepare_manifest.result is not None and prepare_manifest.result.prepared_states:
        return tuple(entry.state_name for entry in prepare_manifest.result.prepared_states)

    state_names = tuple(entry.state_name for entry in prepare_manifest.spec.state_inputs)
    if not state_names:
        raise ValueError("Prepare manifest does not define any state inputs.")
    return state_names


def _resolved_input_record_files(
    spec: AnalyzeWorkflowSpec,
    run_manifest: RunWorkflowManifest,
) -> tuple[StageFileRecord, ...]:
    if spec.input_record_files:
        return tuple(spec.input_record_files)

    if run_manifest.result is None or not run_manifest.result.sample_record_files:
        raise ValueError(
            "Analyze workflow requires 'input_record_files' or a run manifest with sample_record_files."
        )

    json_records = tuple(
        record for record in run_manifest.result.sample_record_files if record.path.suffix.lower() == ".json"
    )
    if not json_records:
        raise ValueError(
            "Run manifest does not contain JSON sample record files suitable for analysis. "
            "Provide 'input_record_files' explicitly."
        )
    return json_records


def _load_window_records(input_records: tuple[StageFileRecord, ...]):
    window_records = []
    for record in input_records:
        window_records.extend(read_fep_window_records(record.path))
    if not window_records:
        raise ValueError("Analyze workflow did not load any FEP window records.")
    return tuple(window_records)


def _load_replica_window_records(input_records: tuple[StageFileRecord, ...]) -> tuple[tuple[object, ...], ...]:
    replica_records = tuple(tuple(read_fep_window_records(record.path)) for record in input_records)
    if not replica_records:
        raise ValueError("Replica analysis requires at least one input record file.")
    for index, records in enumerate(replica_records):
        if not records:
            raise ValueError(f"Replica input record file at index {index} did not contain any FEP window records.")
    return replica_records


def _thermal_energy_kj_mol(spec: AnalyzeWorkflowSpec, run_manifest: RunWorkflowManifest) -> float:
    temperature_kelvin = float(spec.temperature_kelvin or run_manifest.spec.temperature_kelvin or 300.0)
    temperature = temperature_kelvin * unit.kelvin
    return float((unit.MOLAR_GAS_CONSTANT_R * temperature).value_in_unit(unit.kilojoules_per_mole))


def _write_json(path: Path, payload: Mapping[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(_jsonify(payload), indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _jsonify(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if is_dataclass(value):
        # Shallow conversion (not recursive asdict) so _jsonify recurses into each
        # nested dataclass itself -- this lets the per-record None-omission below
        # reach nested FEPSampleRecords. Optional record fields that default to
        # None are omitted when unset, so adding them is purely additive: existing
        # serialized outputs (and golden fixtures) stay byte-identical.
        # FEPSampleRecord.coupling_distance_nm is only populated for the
        # distance-dependent-H12 workflow.
        data = {f.name: getattr(value, f.name) for f in fields(value)}
        if data.get("coupling_distance_nm", False) is None:
            data.pop("coupling_distance_nm", None)
        return _jsonify(data)
    if isinstance(value, Mapping):
        return {str(key): _jsonify(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return [_jsonify(item) for item in value]
    if isinstance(value, list):
        return [_jsonify(item) for item in value]
    return value


def _fmt_optional(value: object) -> str:
    if value is None:
        return ""
    return f"{float(value):.8f}"


__all__ = ["analyze_workflow"]
