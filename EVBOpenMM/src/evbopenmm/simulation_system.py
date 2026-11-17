from __future__ import annotations

from .system import EVBSystem

EVBSimulationSystem = EVBSystem
"""Cluster-side runtime entry point.

``EVBSimulationSystem`` is currently an alias for :class:`EVBSystem`. It exists
as a named home for the "consume a :class:`~evbopenmm.prep.PreparedBundle`,
produce a simulation-ready EVB system" use case — the path that does not need
the OpenFF/RDKit/AmberTools chemistry stack. Reach for it via::

    from evbopenmm import EVBSimulationSystem
    from evbopenmm.prep import PreparedBundle

    bundle = PreparedBundle.from_parameter_bundle_json(path, base_dir=root)
    system = EVBSimulationSystem.from_prepared_bundle(bundle)

The preparation-side counterpart is :class:`evbopenmm.prep.EVBPreparation`.
"""


__all__ = ["EVBSimulationSystem"]
