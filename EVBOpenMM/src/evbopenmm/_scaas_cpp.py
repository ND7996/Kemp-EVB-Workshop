from __future__ import annotations

import ctypes
from pathlib import Path
import os
import platform as sys_platform
import sys
from typing import Tuple

import numpy as np
import openmm


_REFERENCE_UPDATE_SYMBOL = "scaas_reference_update"
_BULK_SET_BONDS_SYMBOL = "scaas_bulk_set_polarization_bonds"
_CUDA_RUNTIME_IS_REGISTERED_SYMBOL = "scaas_cuda_runtime_is_registered"
_CUDA_RUNTIME_DIAGNOSTICS_SYMBOL = "scaas_cuda_runtime_get_diagnostics"
_CUDA_RUNTIME_DIAGNOSTICS_EX_SYMBOL = "scaas_cuda_runtime_get_diagnostics_ex"
_CUDA_RUNTIME_TRIGGER_UPDATE_SYMBOL = "scaas_cuda_runtime_trigger_update"
_CTYPES_LIBRARY: ctypes.CDLL | None = None
_CTYPES_LIBRARY_PATH: Path | None = None
_REFERENCE_UPDATE_FUNCTION: ctypes._CFuncPtr | None = None
_BULK_SET_BONDS_FUNCTION: ctypes._CFuncPtr | None = None
_CUDA_RUNTIME_IS_REGISTERED_FUNCTION: ctypes._CFuncPtr | None = None
_CUDA_RUNTIME_DIAGNOSTICS_FUNCTION: ctypes._CFuncPtr | None = None
_CUDA_RUNTIME_DIAGNOSTICS_EX_FUNCTION: ctypes._CFuncPtr | None = None
_CUDA_RUNTIME_TRIGGER_UPDATE_FUNCTION: ctypes._CFuncPtr | None = None
_CUDA_RUNTIME_DIAGNOSTICS_EX_CHECKED = False


def _platform_library_names() -> Tuple[str, ...]:
    system = sys_platform.system().lower()
    if system == "windows":
        return ("OpenMMSCAAS.dll",)
    if system == "darwin":
        return ("libOpenMMSCAAS.dylib", "OpenMMSCAAS.dylib")
    return ("libOpenMMSCAAS.so", "OpenMMSCAAS.so")


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _dedupe_paths(paths) -> Tuple[Path, ...]:
    seen: set[str] = set()
    ordered: list[Path] = []
    for raw_path in paths:
        path = Path(raw_path)
        key = str(path)
        if key in seen:
            continue
        seen.add(key)
        ordered.append(path)
    return tuple(ordered)


def _active_python_prefix() -> Path:
    return Path(sys.prefix).resolve()


def _openmm_package_root() -> Path:
    return Path(openmm.__file__).resolve().parent


def _default_openmm_plugin_directory() -> Path | None:
    try:
        raw = openmm.Platform.getDefaultPluginsDirectory()
    except openmm.OpenMMException:
        return None
    if not raw:
        return None
    return Path(raw).expanduser().resolve()


def _environment_plugin_directories() -> Tuple[Path, ...]:
    candidates = []
    openmm_plugin_dir = os.environ.get("OPENMM_PLUGIN_DIR")
    if openmm_plugin_dir:
        candidates.append(Path(openmm_plugin_dir))
    default_plugin_dir = _default_openmm_plugin_directory()
    if default_plugin_dir is not None:
        candidates.append(default_plugin_dir)

    package_root = _openmm_package_root()
    site_packages_root = package_root.parent
    prefix = _active_python_prefix()
    system = sys_platform.system().lower()

    candidates.extend(
        [
            site_packages_root / "OpenMM.libs" / "lib" / "plugins",
            site_packages_root / "OpenMM.libs" / "plugins",
            package_root / "plugins",
        ]
    )
    if system == "windows":
        candidates.extend(
            [
                prefix / "Library" / "lib" / "plugins",
                prefix / "Library" / "bin",
            ]
        )
    else:
        candidates.extend(
            [
                prefix / "lib" / "plugins",
                prefix / "plugins",
            ]
        )
    return _dedupe_paths(candidates)


def _environment_library_directories() -> Tuple[Path, ...]:
    candidates = []
    package_root = _openmm_package_root()
    site_packages_root = package_root.parent
    prefix = _active_python_prefix()
    system = sys_platform.system().lower()

    openmm_plugin_dir = os.environ.get("OPENMM_PLUGIN_DIR")
    if openmm_plugin_dir:
        plugin_dir = Path(openmm_plugin_dir)
        candidates.extend([plugin_dir.parent, plugin_dir])
    default_plugin_dir = _default_openmm_plugin_directory()
    if default_plugin_dir is not None:
        candidates.extend([default_plugin_dir.parent, default_plugin_dir])

    candidates.extend(
        [
            site_packages_root / "OpenMM.libs" / "lib",
            site_packages_root / "OpenMM.libs",
        ]
    )
    if system == "windows":
        candidates.extend(
            [
                prefix / "Library" / "lib",
                prefix / "Library" / "bin",
            ]
        )
    else:
        candidates.append(prefix / "lib")
    return _dedupe_paths(candidates)


def _candidate_directories() -> Tuple[Path, ...]:
    custom_dir = os.environ.get("EVBOPENMM_SCAAS_PLUGIN_DIR")
    candidates = []
    if custom_dir:
        candidates.append(Path(custom_dir))

    candidates.extend(_environment_plugin_directories())
    candidates.extend(_environment_library_directories())

    repo = _repo_root()
    candidates.extend(
        [
            repo / "build/scaas-plugin",
            repo / "build/plugins/scaas",
            repo / "plugins/scaas/build",
        ]
    )
    return _dedupe_paths(candidates)


def candidate_plugin_libraries() -> Tuple[Path, ...]:
    """
    Enumerate plausible plugin-library paths for the SCAAS plugin.

    Priority:
    1) ``EVBOPENMM_SCAAS_PLUGIN_LIBRARY`` if set.
    2) Active-environment OpenMM/plugin folders.
    3) Known repository build folders with platform-specific library names.
    """

    explicit = os.environ.get("EVBOPENMM_SCAAS_PLUGIN_LIBRARY")
    if explicit:
        return (Path(explicit),)

    names = _platform_library_names()
    paths: list[Path] = []
    for directory in _candidate_directories():
        for name in names:
            paths.append(directory / name)
    return _dedupe_paths(paths)


def find_plugin_library() -> Path | None:
    """Return first existing candidate plugin library path, if any."""

    for path in candidate_plugin_libraries():
        if path.exists():
            return path
    return None


def load_scaas_plugin(library_path: os.PathLike[str] | str | None = None) -> Path:
    """
    Load the SCAAS plugin library into OpenMM.

    Returns the resolved library path on success.
    """

    if library_path is None:
        resolved = find_plugin_library()
    else:
        resolved = Path(library_path)

    if resolved is None:
        searched = ", ".join(str(p) for p in candidate_plugin_libraries())
        raise FileNotFoundError(f"SCAAS plugin library not found. Searched: {searched}")

    openmm.Platform.loadPluginLibrary(str(resolved))
    return resolved


def _load_ctypes_library(library_path: os.PathLike[str] | str | None = None) -> tuple[ctypes.CDLL, Path]:
    global _CTYPES_LIBRARY
    global _CTYPES_LIBRARY_PATH

    if library_path is None:
        resolved = find_plugin_library()
    else:
        resolved = Path(library_path)

    if resolved is None:
        searched = ", ".join(str(p) for p in candidate_plugin_libraries())
        raise FileNotFoundError(f"SCAAS plugin library not found. Searched: {searched}")

    if _CTYPES_LIBRARY is not None and _CTYPES_LIBRARY_PATH == resolved:
        return _CTYPES_LIBRARY, resolved

    library = ctypes.CDLL(str(resolved))
    _CTYPES_LIBRARY = library
    _CTYPES_LIBRARY_PATH = resolved
    return library, resolved


def probe_scaas_plugin() -> Tuple[bool, str]:
    """
    Probe whether the SCAAS plugin library can be loaded.

    This only checks shared-library loadability.
    """

    path = find_plugin_library()
    if path is None:
        searched = ", ".join(str(p) for p in candidate_plugin_libraries())
        return False, f"SCAAS plugin library not found (searched: {searched})"

    try:
        load_scaas_plugin(path)
        _load_ctypes_library(path)
    except (OSError, openmm.OpenMMException) as exc:
        return False, f"Found '{path}' but failed to load: {exc}"

    return True, f"Loaded plugin library '{path}'."


def probe_reference_updater() -> Tuple[bool, str]:
    """
    Probe whether the native reference updater symbol is available.

    Returns ``(ready, detail)`` where ``ready`` means the symbol exists and the
    library can be loaded.
    """

    ok, detail = probe_scaas_plugin()
    if not ok:
        return False, detail

    library, path = _load_ctypes_library()
    if not hasattr(library, _REFERENCE_UPDATE_SYMBOL):
        return False, f"Loaded '{path}' but symbol '{_REFERENCE_UPDATE_SYMBOL}' is missing."

    return True, f"Loaded '{path}' with symbol '{_REFERENCE_UPDATE_SYMBOL}'."


def probe_bulk_parameter_updater() -> Tuple[bool, str]:
    """
    Probe whether the native bulk bond-parameter updater symbol is available.

    Returns ``(ready, detail)`` where ``ready`` means the symbol exists and the
    library can be loaded.
    """

    ok, detail = probe_scaas_plugin()
    if not ok:
        return False, detail

    library, path = _load_ctypes_library()
    if not hasattr(library, _BULK_SET_BONDS_SYMBOL):
        return False, f"Loaded '{path}' but symbol '{_BULK_SET_BONDS_SYMBOL}' is missing."

    return True, f"Loaded '{path}' with symbol '{_BULK_SET_BONDS_SYMBOL}'."


def probe_cuda_native_runtime() -> Tuple[bool, str]:
    """
    Probe whether native CUDA runtime ownership symbols are available.

    Returns ``(ready, detail)`` where ``ready`` means both symbols exist and
    the plugin can be loaded.
    """

    ok, detail = probe_scaas_plugin()
    if not ok:
        return False, detail

    library, path = _load_ctypes_library()
    missing = []
    for symbol in (_CUDA_RUNTIME_IS_REGISTERED_SYMBOL, _CUDA_RUNTIME_DIAGNOSTICS_SYMBOL):
        if not hasattr(library, symbol):
            missing.append(symbol)
    if missing:
        missing_text = ", ".join(missing)
        return False, f"Loaded '{path}' but missing symbol(s): {missing_text}."

    return True, f"Loaded '{path}' with CUDA native-runtime diagnostics symbols."


def probe_cuda_native_runtime_trigger() -> Tuple[bool, str]:
    """
    Probe whether the native CUDA runtime trigger symbol is available.

    Returns ``(ready, detail)`` where ``ready`` means the plugin can be loaded
    and exposes the direct-update hook used by shared-trajectory replay.
    """

    ok, detail = probe_cuda_native_runtime()
    if not ok:
        return False, detail

    library, path = _load_ctypes_library()
    if not hasattr(library, _CUDA_RUNTIME_TRIGGER_UPDATE_SYMBOL):
        return False, f"Loaded '{path}' but symbol '{_CUDA_RUNTIME_TRIGGER_UPDATE_SYMBOL}' is missing."

    return True, f"Loaded '{path}' with CUDA native-runtime trigger symbol."


def _get_reference_update_function() -> ctypes._CFuncPtr:
    global _REFERENCE_UPDATE_FUNCTION

    if _REFERENCE_UPDATE_FUNCTION is not None:
        return _REFERENCE_UPDATE_FUNCTION

    library, _ = _load_ctypes_library()
    fn = getattr(library, _REFERENCE_UPDATE_SYMBOL)
    fn.argtypes = [
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
    ]
    fn.restype = ctypes.c_int
    _REFERENCE_UPDATE_FUNCTION = fn
    return fn


def _get_bulk_set_bonds_function() -> ctypes._CFuncPtr:
    global _BULK_SET_BONDS_FUNCTION

    if _BULK_SET_BONDS_FUNCTION is not None:
        return _BULK_SET_BONDS_FUNCTION

    library, _ = _load_ctypes_library()
    fn = getattr(library, _BULK_SET_BONDS_SYMBOL)
    fn.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
    ]
    fn.restype = ctypes.c_int
    _BULK_SET_BONDS_FUNCTION = fn
    return fn


def _get_cuda_runtime_is_registered_function() -> ctypes._CFuncPtr:
    global _CUDA_RUNTIME_IS_REGISTERED_FUNCTION

    if _CUDA_RUNTIME_IS_REGISTERED_FUNCTION is not None:
        return _CUDA_RUNTIME_IS_REGISTERED_FUNCTION

    library, _ = _load_ctypes_library()
    fn = getattr(library, _CUDA_RUNTIME_IS_REGISTERED_SYMBOL)
    fn.argtypes = [ctypes.c_void_p]
    fn.restype = ctypes.c_int
    _CUDA_RUNTIME_IS_REGISTERED_FUNCTION = fn
    return fn


def _get_cuda_runtime_diagnostics_function() -> ctypes._CFuncPtr:
    global _CUDA_RUNTIME_DIAGNOSTICS_FUNCTION

    if _CUDA_RUNTIME_DIAGNOSTICS_FUNCTION is not None:
        return _CUDA_RUNTIME_DIAGNOSTICS_FUNCTION

    library, _ = _load_ctypes_library()
    fn = getattr(library, _CUDA_RUNTIME_DIAGNOSTICS_SYMBOL)
    fn.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_longlong),
        ctypes.POINTER(ctypes.c_longlong),
        ctypes.POINTER(ctypes.c_longlong),
        ctypes.POINTER(ctypes.c_longlong),
        ctypes.POINTER(ctypes.c_longlong),
        ctypes.POINTER(ctypes.c_longlong),
        ctypes.POINTER(ctypes.c_longlong),
    ]
    fn.restype = ctypes.c_int
    _CUDA_RUNTIME_DIAGNOSTICS_FUNCTION = fn
    return fn


def _get_cuda_runtime_diagnostics_ex_function() -> ctypes._CFuncPtr | None:
    global _CUDA_RUNTIME_DIAGNOSTICS_EX_FUNCTION
    global _CUDA_RUNTIME_DIAGNOSTICS_EX_CHECKED

    if _CUDA_RUNTIME_DIAGNOSTICS_EX_FUNCTION is not None:
        return _CUDA_RUNTIME_DIAGNOSTICS_EX_FUNCTION
    if _CUDA_RUNTIME_DIAGNOSTICS_EX_CHECKED:
        return None

    library, _ = _load_ctypes_library()
    _CUDA_RUNTIME_DIAGNOSTICS_EX_CHECKED = True
    if not hasattr(library, _CUDA_RUNTIME_DIAGNOSTICS_EX_SYMBOL):
        return None

    fn = getattr(library, _CUDA_RUNTIME_DIAGNOSTICS_EX_SYMBOL)
    fn.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_longlong),
        ctypes.POINTER(ctypes.c_longlong),
        ctypes.POINTER(ctypes.c_longlong),
        ctypes.POINTER(ctypes.c_longlong),
        ctypes.POINTER(ctypes.c_longlong),
        ctypes.POINTER(ctypes.c_longlong),
        ctypes.POINTER(ctypes.c_longlong),
        ctypes.POINTER(ctypes.c_longlong),
        ctypes.POINTER(ctypes.c_longlong),
        ctypes.POINTER(ctypes.c_longlong),
        ctypes.POINTER(ctypes.c_longlong),
    ]
    fn.restype = ctypes.c_int
    _CUDA_RUNTIME_DIAGNOSTICS_EX_FUNCTION = fn
    return fn


def _get_cuda_runtime_trigger_update_function() -> ctypes._CFuncPtr:
    global _CUDA_RUNTIME_TRIGGER_UPDATE_FUNCTION

    if _CUDA_RUNTIME_TRIGGER_UPDATE_FUNCTION is not None:
        return _CUDA_RUNTIME_TRIGGER_UPDATE_FUNCTION

    library, _ = _load_ctypes_library()
    fn = getattr(library, _CUDA_RUNTIME_TRIGGER_UPDATE_SYMBOL)
    fn.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_longlong]
    fn.restype = ctypes.c_int
    _CUDA_RUNTIME_TRIGGER_UPDATE_FUNCTION = fn
    return fn


def _as_int32_array(values: np.ndarray) -> np.ndarray:
    return np.ascontiguousarray(values, dtype=np.int32)


def _as_float64_array(values: np.ndarray) -> np.ndarray:
    return np.ascontiguousarray(values, dtype=np.float64)


def _as_void_pointer_from_swig(obj: object) -> ctypes.c_void_p:
    swig_ptr = getattr(obj, "this", None)
    if swig_ptr is None:
        raise TypeError("OpenMM SWIG object is missing '.this' pointer.")
    return ctypes.c_void_p(int(swig_ptr))


def run_reference_update(
    *,
    positions_nm: np.ndarray,
    oxygen_indices: np.ndarray,
    hydrogen1_indices: np.ndarray,
    hydrogen2_indices: np.ndarray,
    shell_outer_radius_nm: np.ndarray,
    shell_thickness_nm: np.ndarray,
    shell_cstb: np.ndarray,
    center_nm: tuple[float, float, float],
    update_correction: bool,
    charge_correction: bool,
    step: int,
    update_interval: int,
    theta_corr: np.ndarray,
    theta_corr_accumulator: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Run one native reference-update step and return per-molecule parameters."""

    ready, detail = probe_reference_updater()
    if not ready:
        raise RuntimeError(f"Native SCAAS reference updater is unavailable: {detail}")

    positions = _as_float64_array(positions_nm)
    if positions.ndim != 2 or positions.shape[1] != 3:
        raise ValueError("'positions_nm' must have shape (n_atoms, 3).")

    oxygen = _as_int32_array(oxygen_indices)
    h1 = _as_int32_array(hydrogen1_indices)
    h2 = _as_int32_array(hydrogen2_indices)

    n_molecules = int(oxygen.shape[0])
    if int(h1.shape[0]) != n_molecules or int(h2.shape[0]) != n_molecules:
        raise ValueError("Molecule index arrays must all have identical length.")

    shell_outer = _as_float64_array(shell_outer_radius_nm)
    shell_thickness = _as_float64_array(shell_thickness_nm)
    shell_cstb_arr = _as_float64_array(shell_cstb)
    n_shells = int(shell_outer.shape[0])
    if int(shell_thickness.shape[0]) != n_shells or int(shell_cstb_arr.shape[0]) != n_shells:
        raise ValueError("Shell arrays must all have identical length.")

    theta_corr_arr = _as_float64_array(theta_corr)
    theta_corr_acc_arr = _as_float64_array(theta_corr_accumulator)
    if int(theta_corr_arr.shape[0]) != n_shells or int(theta_corr_acc_arr.shape[0]) != n_shells:
        raise ValueError("'theta_corr' arrays must match shell count.")

    out_active = np.zeros(n_molecules, dtype=np.float64)
    out_targets = np.zeros(n_molecules, dtype=np.float64)
    out_corrections = np.zeros(n_molecules, dtype=np.float64)

    fn = _get_reference_update_function()
    status = int(
        fn(
            ctypes.c_int(n_molecules),
            oxygen.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            h1.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            h2.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            ctypes.c_int(n_shells),
            shell_outer.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            shell_thickness.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            shell_cstb_arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ctypes.c_int(int(positions.shape[0])),
            positions.ravel().ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ctypes.c_double(float(center_nm[0])),
            ctypes.c_double(float(center_nm[1])),
            ctypes.c_double(float(center_nm[2])),
            ctypes.c_int(1 if update_correction else 0),
            ctypes.c_int(1 if charge_correction else 0),
            ctypes.c_int(int(step)),
            ctypes.c_int(int(update_interval)),
            theta_corr_arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            theta_corr_acc_arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            out_active.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            out_targets.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            out_corrections.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        )
    )
    if status != 0:
        raise RuntimeError(f"Native SCAAS reference updater failed with status code {status}.")

    if theta_corr_arr is not theta_corr:
        theta_corr[...] = theta_corr_arr
    if theta_corr_acc_arr is not theta_corr_accumulator:
        theta_corr_accumulator[...] = theta_corr_acc_arr

    return out_active, out_targets, out_corrections


def run_bulk_set_polarization_bonds(
    *,
    force: openmm.CustomCompoundBondForce,
    oxygen_indices: np.ndarray,
    hydrogen1_indices: np.ndarray,
    hydrogen2_indices: np.ndarray,
    active: np.ndarray,
    targets: np.ndarray,
    corrections: np.ndarray,
) -> None:
    """Apply per-bond polarization parameters through one native bulk call."""

    ready, detail = probe_bulk_parameter_updater()
    if not ready:
        raise RuntimeError(f"Native SCAAS bulk updater is unavailable: {detail}")

    oxygen = _as_int32_array(oxygen_indices)
    h1 = _as_int32_array(hydrogen1_indices)
    h2 = _as_int32_array(hydrogen2_indices)
    active_arr = _as_float64_array(active)
    targets_arr = _as_float64_array(targets)
    corrections_arr = _as_float64_array(corrections)

    n_bonds = int(oxygen.shape[0])
    if int(h1.shape[0]) != n_bonds or int(h2.shape[0]) != n_bonds:
        raise ValueError("Bond particle-index arrays must all have identical length.")
    if int(active_arr.shape[0]) != n_bonds or int(targets_arr.shape[0]) != n_bonds or int(corrections_arr.shape[0]) != n_bonds:
        raise ValueError("Bond parameter arrays must all match particle-index array length.")
    if int(force.getNumBonds()) != n_bonds:
        raise ValueError("'force' bond count does not match provided array length.")

    fn = _get_bulk_set_bonds_function()
    status = int(
        fn(
            _as_void_pointer_from_swig(force),
            ctypes.c_int(n_bonds),
            oxygen.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            h1.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            h2.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            active_arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            targets_arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            corrections_arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        )
    )
    if status != 0:
        raise RuntimeError(f"Native SCAAS bulk updater failed with status code {status}.")


def is_cuda_native_runtime_registered(*, force: openmm.CustomCompoundBondForce) -> bool:
    """Return whether CUDA native-runtime ownership is active for this force."""

    ready, _ = probe_cuda_native_runtime()
    if not ready:
        return False

    fn = _get_cuda_runtime_is_registered_function()
    status = int(fn(_as_void_pointer_from_swig(force)))
    return bool(status == 1)


def _force_global_parameter_defaults(force: object) -> dict[str, float]:
    try:
        n_parameters = int(force.getNumGlobalParameters())
    except (AttributeError, TypeError, openmm.OpenMMException):
        return {}
    payload: dict[str, float] = {}
    for index in range(n_parameters):
        try:
            name = str(force.getGlobalParameterName(index))
            value = float(force.getGlobalParameterDefaultValue(index))
        except (AttributeError, TypeError, openmm.OpenMMException):
            continue
        payload[name] = value
    return payload


def _derive_charge_correction_diagnostics(
    *,
    force: object,
    report_count: int,
    last_step: int | None,
) -> dict[str, object]:
    global_parameters = _force_global_parameter_defaults(force)
    raw_update_interval = global_parameters.get("scaas_update_interval")
    update_interval = None
    if raw_update_interval is not None:
        update_interval = max(1, int(round(float(raw_update_interval))))
    charge_correction_enabled = bool(global_parameters.get("scaas_charge_correction", 0.0) > 0.5)

    if (
        not charge_correction_enabled
        or update_interval is None
        or last_step is None
        or last_step <= 0
        or report_count <= 0
    ):
        return {
            "charge_correction_enabled": charge_correction_enabled,
            "update_interval_steps": update_interval,
            "charge_correction_updates": 0,
            "last_charge_correction_update_step": None,
        }

    correction_count = max(0, int(last_step) // int(update_interval))
    correction_count = min(correction_count, int(report_count))
    if correction_count <= 0:
        last_correction_step = None
    else:
        last_correction_step = correction_count * int(update_interval)
    return {
        "charge_correction_enabled": charge_correction_enabled,
        "update_interval_steps": update_interval,
        "charge_correction_updates": correction_count,
        "last_charge_correction_update_step": last_correction_step,
    }


def get_cuda_native_runtime_diagnostics(
    *,
    force: openmm.CustomCompoundBondForce,
    max_shells: int = 256,
) -> dict[str, object]:
    """Fetch native CUDA runtime diagnostics for a tagged SCAAS force."""

    ready, detail = probe_cuda_native_runtime()
    if not ready:
        raise RuntimeError(f"Native CUDA runtime diagnostics are unavailable: {detail}")
    if max_shells <= 0:
        raise ValueError("'max_shells' must be > 0.")

    fn = _get_cuda_runtime_diagnostics_function()
    force_ptr = _as_void_pointer_from_swig(force)

    shell_buffer = int(max_shells)
    while True:
        theta_corr = np.zeros(shell_buffer, dtype=np.float64)
        shell_count = ctypes.c_int(0)
        report_calls = ctypes.c_longlong(0)
        prime_calls = ctypes.c_longlong(0)
        last_update_step = ctypes.c_longlong(-1)
        compute_ns_total = ctypes.c_longlong(0)
        set_params_ns_total = ctypes.c_longlong(0)
        total_update_ns = ctypes.c_longlong(0)
        execute_calls = ctypes.c_longlong(0)
        globals_upload_ns_total = ctypes.c_longlong(0)
        posq_download_ns_total = ctypes.c_longlong(0)
        posq_unpack_ns_total = ctypes.c_longlong(0)
        param_pack_ns_total = ctypes.c_longlong(0)

        fn_ex = _get_cuda_runtime_diagnostics_ex_function()
        if fn_ex is None:
            status = int(
                fn(
                    force_ptr,
                    ctypes.c_int(shell_buffer),
                    theta_corr.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                    ctypes.byref(shell_count),
                    ctypes.byref(report_calls),
                    ctypes.byref(prime_calls),
                    ctypes.byref(last_update_step),
                    ctypes.byref(compute_ns_total),
                    ctypes.byref(set_params_ns_total),
                    ctypes.byref(total_update_ns),
                    ctypes.byref(execute_calls),
                )
            )
        else:
            status = int(
                fn_ex(
                    force_ptr,
                    ctypes.c_int(shell_buffer),
                    theta_corr.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                    ctypes.byref(shell_count),
                    ctypes.byref(report_calls),
                    ctypes.byref(prime_calls),
                    ctypes.byref(last_update_step),
                    ctypes.byref(compute_ns_total),
                    ctypes.byref(set_params_ns_total),
                    ctypes.byref(total_update_ns),
                    ctypes.byref(execute_calls),
                    ctypes.byref(globals_upload_ns_total),
                    ctypes.byref(posq_download_ns_total),
                    ctypes.byref(posq_unpack_ns_total),
                    ctypes.byref(param_pack_ns_total),
                )
            )

        if status == 2:
            return {}
        if status != 0:
            raise RuntimeError(f"Native CUDA runtime diagnostics failed with status code {status}.")

        n_shells = int(shell_count.value)
        if n_shells <= shell_buffer:
            break
        shell_buffer = n_shells

    report_count = int(report_calls.value)
    prime_count = int(prime_calls.value)
    update_divisor = float(report_count) if report_count > 0 else 1.0
    last_step_value = int(last_update_step.value)
    last_step = None if last_step_value < 0 else last_step_value
    theta_values = tuple(float(value) for value in theta_corr[:n_shells].tolist())

    profiling = {
        "compute_stage": "cpp_cuda_native_runtime",
        "report_calls": report_count,
        "prime_calls": prime_count,
        "execute_calls": int(execute_calls.value),
        "compute_ns_total": int(compute_ns_total.value),
        "set_bond_parameters_ns_total": int(set_params_ns_total.value),
        "globals_upload_ns_total": int(globals_upload_ns_total.value),
        "posq_download_ns_total": int(posq_download_ns_total.value),
        "posq_unpack_ns_total": int(posq_unpack_ns_total.value),
        "param_pack_ns_total": int(param_pack_ns_total.value),
        "update_parameters_in_context_ns_total": 0,
        "total_report_ns": int(total_update_ns.value),
        "compute_ns_mean": float(compute_ns_total.value) / update_divisor if report_count > 0 else 0.0,
        "set_bond_parameters_ns_mean": float(set_params_ns_total.value) / update_divisor if report_count > 0 else 0.0,
        "globals_upload_ns_mean": float(globals_upload_ns_total.value) / update_divisor if report_count > 0 else 0.0,
        "posq_download_ns_mean": float(posq_download_ns_total.value) / update_divisor if report_count > 0 else 0.0,
        "posq_unpack_ns_mean": float(posq_unpack_ns_total.value) / update_divisor if report_count > 0 else 0.0,
        "param_pack_ns_mean": float(param_pack_ns_total.value) / update_divisor if report_count > 0 else 0.0,
        "update_parameters_in_context_ns_mean": 0.0,
        "total_report_ns_mean": float(total_update_ns.value) / update_divisor if report_count > 0 else 0.0,
    }
    payload = {
        "polarisation_updates": report_count,
        "polarisation_runtime_updates": report_count,
        "last_polarisation_update_step": last_step,
        "theta_corr": theta_values,
        "profiling": profiling,
    }
    payload.update(
        _derive_charge_correction_diagnostics(
            force=force,
            report_count=report_count,
            last_step=last_step,
        )
    )
    return payload


def trigger_cuda_native_runtime_update(
    *,
    force: openmm.CustomCompoundBondForce,
    context: openmm.Context,
    step: int,
) -> None:
    """Trigger one native CUDA SCAAS update on the live runtime kernel."""

    ready, detail = probe_cuda_native_runtime_trigger()
    if not ready:
        raise RuntimeError(f"Native CUDA runtime trigger is unavailable: {detail}")

    fn = _get_cuda_runtime_trigger_update_function()
    status = int(
        fn(
            _as_void_pointer_from_swig(force),
            _as_void_pointer_from_swig(context),
            ctypes.c_longlong(int(step)),
        )
    )
    if status == 0:
        return
    if status == 2:
        raise RuntimeError("Native CUDA runtime trigger is not registered for the provided force/context.")
    raise RuntimeError(f"Native CUDA runtime trigger failed with status code {status}.")
