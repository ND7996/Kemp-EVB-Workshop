"""
EVB OpenMM - input reader

Public API for loading one or more EVB PDB states into OpenMM objects.

Example
-------
>>> from evbopenmm.input_reader import read_evb_pdbs
>>> state_a, state_b = read_evb_pdbs("state_a.pdb", "state_b.pdb")
>>> top_a, top_b = state_a.topology, state_b.topology
>>> pos_a, pos_b = state_a.positions, state_b.positions
"""

from .system_reader import EVBState, read_evb_pdbs, read_pdb_state

__all__ = ["EVBState", "read_evb_pdbs", "read_pdb_state"]
