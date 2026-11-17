from __future__ import annotations

import argparse
from pathlib import Path
import sys

from .execute import execute_analyze, execute_manifest, execute_prepare, execute_project, execute_run
from .executor import execute_manifest_with_resume
from .models import AnalyzeWorkflowManifest, PrepareWorkflowManifest, ProjectWorkflowManifest, RunWorkflowManifest
from .project import (
    PROJECT_RESTART_POLICY_CASCADE,
    VALID_PROJECT_RECOVERY_START_STAGES,
    VALID_PROJECT_RESTART_POLICIES,
)
from .project import recover_project_manifest
from .status import collect_workflow_status


def build_parser() -> argparse.ArgumentParser:
    """Build the minimal manifest-first workflow CLI parser."""

    parser = argparse.ArgumentParser(
        prog="evbopenmm-workflow",
        description=(
            "Execute EVBOpenMM workflow manifests. Primary commands are "
            "prepare, run, analyze, project, recover, and status; execute and "
            "resume are advanced dispatch/recovery helpers."
        ),
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    for command_name, help_text in (
        ("prepare", "Primary: execute a prepare manifest."),
        ("run", "Primary: execute a run manifest."),
        ("analyze", "Primary: execute an analyze manifest."),
        ("project", "Primary: execute a chained prepare -> run -> analyze project manifest."),
        ("recover", "Primary: recover a chained project manifest from status metadata."),
        ("status", "Primary: print workflow status for a stage or project manifest."),
        ("execute", "Advanced: dispatch any supported manifest type."),
        ("resume", "Advanced: rerun one stage manifest with stage-local restart handling."),
    ):
        subparser = subparsers.add_parser(command_name, help=help_text)
        subparser.add_argument("manifest", type=Path, help="Path to the workflow manifest JSON file.")
        if command_name in {"project", "resume"}:
            subparser.add_argument(
                "--force-rerun",
                action="store_true",
                help="Ignore any existing stage-status file and rerun the manifest.",
            )
        if command_name == "project":
            subparser.add_argument(
                "--restart-policy",
                choices=tuple(sorted(VALID_PROJECT_RESTART_POLICIES)),
                default=PROJECT_RESTART_POLICY_CASCADE,
                help="Project-level restart policy for downstream stage reuse.",
            )
        if command_name == "recover":
            recovery_group = subparser.add_mutually_exclusive_group()
            recovery_group.add_argument(
                "--from-stage",
                choices=tuple(sorted(VALID_PROJECT_RECOVERY_START_STAGES)),
                default=None,
                help="Force recovery to rerun from the selected stage onward.",
            )
            recovery_group.add_argument(
                "--analyze-label",
                default=None,
                help="Force recovery to rerun only the named analyze stage label.",
            )

    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the minimal manifest-first EVB workflow CLI."""

    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "status":
            print(collect_workflow_status(args.manifest).to_json())
            return 0

        output_manifest = _execute_command(
            args.command,
            args.manifest,
            force_rerun=bool(getattr(args, "force_rerun", False)),
            restart_policy=str(getattr(args, "restart_policy", PROJECT_RESTART_POLICY_CASCADE)),
            start_stage=getattr(args, "from_stage", None),
            analyze_label=getattr(args, "analyze_label", None),
        )
    except Exception as exc:  # noqa: BLE001 - CLI entry point: report and return non-zero
        print(f"error: {exc}", file=sys.stderr)
        return 1

    output_path = _output_manifest_path(output_manifest)
    print(output_path)
    return 0


def _execute_command(
    command: str,
    manifest_path: Path,
    *,
    force_rerun: bool = False,
    restart_policy: str = PROJECT_RESTART_POLICY_CASCADE,
    start_stage: str | None = None,
    analyze_label: str | None = None,
) -> PrepareWorkflowManifest | RunWorkflowManifest | AnalyzeWorkflowManifest | ProjectWorkflowManifest:
    if command == "prepare":
        return execute_prepare(manifest_path)
    if command == "run":
        return execute_run(manifest_path)
    if command == "analyze":
        return execute_analyze(manifest_path)
    if command == "project":
        return execute_project(manifest_path, force_rerun=force_rerun, restart_policy=restart_policy)
    if command == "recover":
        output_manifest, _ = recover_project_manifest(
            manifest_path,
            start_stage=start_stage,
            analyze_label=analyze_label,
        )
        return output_manifest
    if command == "execute":
        return execute_manifest(manifest_path)
    if command == "resume":
        output_manifest, _ = execute_manifest_with_resume(manifest_path, force_rerun=force_rerun)
        return output_manifest
    raise ValueError(f"Unsupported command {command!r}.")


def _output_manifest_path(
    manifest: PrepareWorkflowManifest | RunWorkflowManifest | AnalyzeWorkflowManifest | ProjectWorkflowManifest,
) -> Path:
    if manifest.result is None or manifest.result.manifest_path is None:
        raise ValueError("Workflow execution did not produce an output manifest path.")
    return manifest.result.manifest_path


__all__ = ["build_parser", "main"]
