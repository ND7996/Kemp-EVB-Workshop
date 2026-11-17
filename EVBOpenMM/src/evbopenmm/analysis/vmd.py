"""Generate a self-contained VMD session (TCL) to inspect an EVB trajectory.

An EVB mapping-FEP / MD run writes a DCD plus a topology PDB.  A plain VMD load is
unsatisfying for EVB because (1) VMD guesses bonds from the *first* frame, so the
breaking bond stays drawn and the forming bond never appears, and (2) the reaction
progress is invisible.  :func:`write_vmd_session` emits a ``.tcl`` you open with
``vmd -e session.tcl`` that solves both, with every EVB-specific feature behind its
own keyword so callers can turn each on/off (project defaults are decided elsewhere):

* ``dynamic_bonds``     -- draw the reactive atoms with a ``DynamicBonds`` rep, so the
                           forming/breaking bonds appear and vanish *as the geometry
                           changes* along the trajectory (the key EVB view).
* ``color_by_reaction`` -- measure the forming-bond distance per frame in TCL and store
                           it in the ``user`` field, colouring the reactive atoms by it
                           -- a self-contained reaction-progress colour, no sidecar data.
* ``all_waters``        -- show every water in the droplet (default; the SCAAS sphere at a glance).
* ``water_shell``       -- show (and auto-update) waters within a cutoff of the core.
* ``show_ions``         -- CPK non-protein ions (e.g. catalytic Ca2+).
* ``show_sphere``       -- outline the SCAAS droplet boundary (caller supplies centre/radius).

The function is pure string generation (no OpenMM/VMD import), so it is trivially
testable and usable offline from saved run outputs.  Reactive atoms are selected BY
INDEX -- droplet builds renumber residues, so name/resid selections are unreliable.
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Iterable, Mapping, Sequence

__all__ = ["write_vmd_session"]


def _as_index_list(reactive_atom_indices) -> list[int]:
    if reactive_atom_indices is None:
        return []
    if isinstance(reactive_atom_indices, Mapping):
        return [int(v) for v in reactive_atom_indices.values()]
    return [int(i) for i in reactive_atom_indices]


def _as_pairs(forming_breaking_pairs) -> list[tuple[int, int]]:
    if not forming_breaking_pairs:
        return []
    pairs = []
    for p in forming_breaking_pairs:
        i, j = p
        pairs.append((int(i), int(j)))
    return pairs


def write_vmd_session(
    output_tcl: str | os.PathLike,
    topology_file: str | os.PathLike,
    dcd_file: str | os.PathLike,
    *,
    reactive_atom_indices: Sequence[int] | Mapping[str, int] | None = None,
    forming_breaking_pairs: Iterable[tuple[int, int]] | None = None,
    dynamic_bonds: bool = True,
    dynamic_bond_cutoff_a: float = 2.2,
    color_by_reaction: bool = True,
    reaction_pair: tuple[int, int] | None = None,
    protein_cartoon: bool = True,
    all_waters: bool = True,
    water_shell: bool = True,
    water_within_a: float = 6.0,
    show_ions: bool = True,
    show_sphere: bool = False,
    sphere_center_nm: Sequence[float] | None = None,
    sphere_radius_a: float | None = None,
    relative_paths: bool = True,
    title: str = "EVB trajectory",
) -> str:
    """Write a VMD TCL session for an EVB trajectory and return the TCL text.

    Parameters
    ----------
    output_tcl, topology_file, dcd_file
        Destination ``.tcl`` and the run's topology PDB + DCD.  When ``relative_paths``
        the loader refers to the trajectory by basename, so the session is portable if
        opened from the run directory.
    reactive_atom_indices
        The reactive-core atoms to highlight, as a list of 0-based topology indices or a
        ``{name: index}`` map.  Selected by index (numbering-proof under droplet renumber).
    forming_breaking_pairs
        ``(i, j)`` index pairs of the forming/breaking bonds.  Drawn via ``DynamicBonds``
        (see ``dynamic_bonds``); the first pair is the default reaction coordinate.
    dynamic_bonds, dynamic_bond_cutoff_a
        Draw a ``DynamicBonds`` rep (cutoff in A) over the reactive atoms so bonds appear
        / disappear as the reaction proceeds.
    color_by_reaction, reaction_pair
        Colour the reactive atoms by the per-frame distance of ``reaction_pair`` (defaults
        to the first forming/breaking pair), computed in TCL into the ``user`` field.
    all_waters
        Draw every water in the system as thin lines.  On by default: in a SCAAS droplet this is
        what shows the sphere holding together (or evaporating) and whether the radius really
        covers the reacting complex.  Independent of ``water_shell`` -- both reps can coexist and
        either can be toggled off in VMD; the session prints the ``mol showrep`` line for each.
    protein_cartoon, water_shell, water_within_a, show_ions
        Optional extra representations, each independently toggled.
    show_sphere, sphere_center_nm, sphere_radius_a
        Outline the SCAAS droplet boundary; centre in nm, radius in A.  Both required when
        ``show_sphere``.
    relative_paths, title
        Portability and the molecule label.

    Returns
    -------
    str
        The TCL text (also written to ``output_tcl``).
    """
    ridx = _as_index_list(reactive_atom_indices)
    pairs = _as_pairs(forming_breaking_pairs)
    if any(i < 0 for i in ridx) or any(i < 0 or j < 0 for i, j in pairs):
        raise ValueError("reactive/pair indices must be non-negative topology indices")
    if dynamic_bond_cutoff_a <= 0 or water_within_a <= 0:
        raise ValueError("cutoff distances must be positive")
    if show_sphere and (sphere_center_nm is None or sphere_radius_a is None):
        raise ValueError("show_sphere requires both sphere_center_nm and sphere_radius_a")

    core = sorted(set(ridx) | {a for p in pairs for a in p})
    core_sel = " ".join(str(i) for i in core)

    if reaction_pair is None and pairs:
        reaction_pair = pairs[0]

    def _path(p):
        return os.path.basename(str(p)) if relative_paths else str(p)

    L: list[str] = []
    L.append(f"# VMD session for an EVB trajectory -- generated by evbopenmm.analysis.vmd")
    L.append(f"# {title}")
    if core:
        L.append(f"# reactive atoms (0-based index): {core_sel}")
    for a, b in pairs:
        L.append(f"# forming/breaking pair: {a} {b}")
    L.append(f"mol new {_path(topology_file)} type pdb waitfor all")
    L.append(f"mol addfile {_path(dcd_file)} type dcd waitfor all")
    L.append('mol rename top "%s"' % title)
    L.append("mol delrep 0 top")

    rep = 0

    def addrep(representation: str, selection: str, *, color: str | None = None,
               selupdate: bool = False):
        nonlocal rep
        L.append(f"mol representation {representation}")
        L.append("mol color %s" % (color or "Name"))
        L.append("mol selection {%s}" % selection)
        L.append("mol addrep top")
        if selupdate:
            L.append(f"mol selupdate {rep} top on")
        this = rep
        rep += 1
        return this

    water_rep = shell_rep = None

    if protein_cartoon:
        addrep("NewCartoon", "protein", color="Structure")

    if core:
        addrep("Licorice 0.3 12", f"same residue as index {core_sel}")
        core_rep = addrep("VDW 0.4 20", f"index {core_sel}")
    else:
        core_rep = None

    if dynamic_bonds and core:
        addrep(f"DynamicBonds {dynamic_bond_cutoff_a:.2f} 0.15 12", f"index {core_sel}")

    if show_ions:
        addrep("CPK 1.0 0.0 20", "ions or (resname CA Ca CAL and not protein)")

    # Two independent water reps so both renderings are available in one session and either can be
    # switched off from the GUI (or with the `mol showrep` lines printed at the end): the whole
    # droplet -- which is what you want in a SCAAS sphere, to see the boundary hold and judge whether
    # the radius actually covers the reacting complex -- and the close shell that follows the core.
    if all_waters:
        water_rep = addrep("Lines 1", "water")
    if water_shell and core:
        shell_rep = addrep("Licorice 0.1 12 12",
                           f"water and same residue as (within {water_within_a:.1f} of index {core_sel})",
                           selupdate=True)

    # reaction-progress colouring: measure the forming-bond distance each frame -> user field.
    if color_by_reaction and reaction_pair is not None and core_rep is not None:
        i, j = reaction_pair
        L.append("")
        L.append("# colour the reactive atoms by the per-frame forming-bond distance (reaction progress)")
        L.append(f"set _rp_i {i}")
        L.append(f"set _rp_j {j}")
        L.append('set _core [atomselect top "index %s"]' % core_sel)
        L.append("set _nf [molinfo top get numframes]")
        L.append("set _dmin 1e9; set _dmax -1e9")
        L.append("for {set f 0} {$f < $_nf} {incr f} {")
        L.append("    set d [measure bond [list $_rp_i $_rp_j] frame $f]")
        L.append("    $_core frame $f")
        L.append("    $_core set user $d")
        L.append("    if {$d < $_dmin} {set _dmin $d}")
        L.append("    if {$d > $_dmax} {set _dmax $d}")
        L.append("}")
        L.append(f"mol modcolor {core_rep} top User")
        L.append(f"mol scaleminmax top {core_rep} $_dmin $_dmax")
        L.append("$_core delete")

    if show_sphere:
        cx, cy, cz = (float(c) * 10.0 for c in sphere_center_nm)  # nm -> A
        L.append("")
        L.append("# SCAAS droplet boundary")
        L.append("draw color blue")
        L.append("draw material Transparent")
        L.append(f"draw sphere {{{cx:.3f} {cy:.3f} {cz:.3f}}} radius {float(sphere_radius_a):.2f} resolution 30")

    if water_rep is not None or shell_rep is not None:
        L.append("")
        L.append("# Two water views are loaded; switch either off without hunting for rep numbers.")
        if water_rep is not None:
            L.append(f"# whole droplet   : mol showrep top {water_rep} off   (on: ... on)")
        if shell_rep is not None:
            L.append(f"# core hydration  : mol showrep top {shell_rep} off   (on: ... on)")

    L.append("")
    L.append("display projection Orthographic")
    L.append("axes location Off")
    L.append("color Display Background white")
    L.append("mol showrep top 0 on")
    L.append("animate goto 0")

    text = "\n".join(L) + "\n"
    Path(output_tcl).write_text(text)
    return text
