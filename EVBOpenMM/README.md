# EVBOpenMM: Empirical Valence Bond (EVB) Implementation in OpenMM

## Overview

The **Empirical Valence Bond (EVB)** method provides a physically transparent and computationally efficient way to model **chemical reactions** within molecular dynamics simulations.  
Instead of introducing new quantum–mechanical terms, EVB represents a reacting system as a **superposition of classical force fields** (“diabatic states”) corresponding to different bonding topologies—typically a **reactant** and a **product** state.

At every molecular geometry, the method constructs a small **Hamiltonian matrix** that couples these diabatic states through an empirical interaction term. Diagonalizing this Hamiltonian yields an **adiabatic potential energy surface (PES)** that smoothly connects reactants and products and can be used for molecular dynamics (MD), free-energy sampling, or transition-state analysis.

## Requirements

- Python ≥ 3.8
- OpenMM ≥ 8.2
- NumPy

Optional extras: `scipy` (for fitting), `openff-toolkit` / `openmmforcefields` / `rdkit` (for small-molecule parameterization). See `pyproject.toml` `[project.optional-dependencies]`.

## Workflow CLI

The manifest-driven workflow layer introduced in Phase B is documented in
[docs/workflow_user_guide.md](docs/workflow_user_guide.md).

For cluster-oriented usage, the intended pattern is:

1. run `prepare` on a machine with the chemistry-stack dependencies available
2. copy the prepared bundle (`prepared.pdb`, `state_system.xml`,
   `parameter_bundle.json`, manifests)
3. run `run` and `analyze` on the target machine with the lighter OpenMM +
   EVBOpenMM runtime stack

Prepared bundles are a frozen force-field contract. If the pure-state
`createSystem(...)` setup changes, re-run `prepare` rather than editing the
runtime stage in place.

For the main user-facing path, prefer:

- `prepare`
- `run`
- `analyze`
- `project`
- `status`
- `recover`

The generic `execute` and low-level `resume` commands remain available as
advanced helpers.

## Tutorials

Runnable tutorial material lives under [tutorials/](tutorials/).

The current notebook-based tutorial is:

- `tutorials/scaas_mapping_fep/ch3cl_cl_evb_scaas_mapping_fep_walkthrough/ch3cl_cl_evb_scaas_mapping_fep_walkthrough.ipynb`

---

## Theory

### 1. The EVB Hamiltonian

Consider a reaction that can be described by two valence-bond states:  
state **1** (reactant) and state **2** (product).  
At any molecular configuration \( x \), the EVB Hamiltonian is defined as

$$
\mathbf{H}(x) =
\begin{pmatrix}
V_1(x) + \Delta_1 & H_{12} \\
H_{12} & V_2(x) + \Delta_2
\end{pmatrix}
$$

where

- $V_1(x)$ and $V_2(x)$ are the **molecular-mechanics (MM) potential energies** for states 1 and 2 at geometry $x$;  
  these energies are obtained from standard force fields with different bonded parameters, charges, or Lennard-Jones terms reflecting the bonding pattern of each state.  
- $\Delta_1$ and $\Delta_2$ are **constant energy offsets** that shift the relative stability of the two states and reproduce the desired reaction free energy ($\Delta_2 - \Delta_1 = \Delta G^\circ$).  
- $H_{12}$ is the **coupling constant**, an empirical parameter controlling the strength of the mixing between states.  
  A larger $H_{12}$ corresponds to stronger coupling and therefore a **lower reaction barrier**.

---

### 2. Adiabatic Energies

The adiabatic energies are obtained by **diagonalizing** the Hamiltonian:

$$
E_{\pm}(x) = \frac{1}{2} \Big[ (V_1 + \Delta_1) + (V_2 + \Delta_2) \pm \sqrt{ \big[(V_1 + \Delta_1) - (V_2 + \Delta_2)\big]^2 + 4H_{12}^2 } \Big] $$

- The **lower branch** $E_{-}(x)$ represents the **ground-state EVB surface** and is used as the potential energy for dynamics.  
- The **upper branch** $E_{+}(x)$ corresponds to the first excited EVB surface, seldom used in chemistry but occasionally informative for electron-transfer or proton-transfer coupling analyses.

Because the diagonalization is analytic, both the **energy and its gradient** are continuous functions of $x$, enabling stable MD integration.

---

### 3. The Energy-Gap Reaction Coordinate

A convenient coordinate describing the progress of the reaction is the **energy gap** between diabatic states:

$$
\Delta V(x) = (V_2 + \Delta_2) - (V_1 + \Delta_1)
$$

- When $\Delta V(x) \ll 0$, the system is reactant-like.  
- When $\Delta V(x) \gg 0$, it is product-like.  
- The **transition region** occurs near $\Delta V(x) \approx 0$, where the two states have similar energy and the coupling $H_{12}$ induces strong mixing.

This coordinate is particularly useful for **umbrella sampling** or **free-energy perturbation** calculations, as it provides a scalar reaction coordinate even for complex molecular geometries.

---

### 4. Diabatic Populations

At each configuration, the instantaneous contribution of each diabatic state to the adiabatic ground state can be computed as

$$
c_1^2 =
\frac{1}{2}
\left(
1 + \frac{\Delta V}{\sqrt{\Delta V^2 + 4H_{12}^2}}
\right),
\qquad
c_2^2 = 1 - c_1^2
$$

These quantities ($c_1^2$, $c_2^2$) are the **fractional populations** of the reactant and product diabatic states, respectively.  
They can be monitored along an MD trajectory to visualize how the system transitions between the two bonding patterns.

---

### 5. Physical Meaning of the Parameters

| Parameter | Role | Typical Determination |
|------------|------|-----------------------|
| $V_1(x), V_2(x)$ | Diabatic MM energies (reactant and product) | Calculated directly from classical force fields |
| $\Delta_1, \Delta_2$ | Energy shifts (control thermodynamics) | Fitted to reproduce $\Delta G^\circ$ or QM reference |
| $H_{12}$ | Coupling constant (controls barrier height) | Tuned to match QM barrier or kinetic data |

In practice, $\Delta_i$ and $H_{12}$ are adjusted so that the EVB free-energy profile reproduces **quantum-mechanical reaction energetics** or **experimental activation parameters**.

---

## Usage (Gas Phase, Two-State EVB)

### 1. Run Mapping FEP and Build EVB Profile

```python
from evbopenmm import EVBSystem

evb = EVBSystem(
    reactant="state_a.pdb",
    product="state_b.pdb",
    force_field="amber14",
)

analysis = evb.run_mapping_fep_evb_profile(
    lambda_start=1.0,
    lambda_end=0.0,
    lambda_step=-0.02,
    equilibration_steps=2000,
    production_steps=10000,
    sample_interval=10,
    n_gap_bins=60,
    min_points_per_bin=20,
    alpha_kcal_mol=0.0,
    hij_kcal_mol=5.0,
    output_directory="results/single_run",
)

print(analysis["thermochemistry"]["dg_dagger_kcal_mol"])
print(analysis["thermochemistry"]["dg0_kcal_mol"])
```

`run_mapping_fep_evb_profile(...)` writes:
- `fep_windows.csv`
- `fep_edges.csv`
- `umbrella_profile.csv`
- `evb_profile.csv`
- `evb_thermochemistry.txt`

### 1b. Optional: Solute Positional Restraints During Pre-equilibration

You can pre-equilibrate with harmonic positional restraints on solute atoms and
release them across equilibration stages (default schedule in
`equilibrate_schedule`: 100 steps at 0.5 fs/25 K, then 100-step windows at
1.0 fs with 50..300 K).

```python
analysis = evb.run_mapping_fep_umbrella(
    lambda_start=0.99,
    lambda_end=0.01,
    lambda_step=-0.01,
    production_steps=500,
    sample_interval=10,
    pre_equilibration=True,
    pre_equilibration_solute_restraints=True,
    pre_equilibration_solute_restraint_k_start=2000.0,  # kJ/mol/nm^2
    pre_equilibration_solute_restraint_k_end=0.0,
)
```

Notes:
- By default, restrained atoms are auto-selected as non-water heavy atoms.
- Override selection with `pre_equilibration_solute_restraint_atom_indices`.
- Override the release profile with `pre_equilibration_solute_restraint_k_schedule`.

### 2. Aggregate Multiple Replica Profiles

```python
from evbopenmm import EVBSystem

evb = EVBSystem(reactant="state_a.pdb", product="state_b.pdb")

replica_window_records = []
for seed in (1001, 1002, 1003):
    records = evb.collect_mapping_fep_samples(
        lambda_start=1.0,
        lambda_end=0.0,
        lambda_step=-0.02,
        equilibration_steps=2000,
        production_steps=10000,
        sample_interval=10,
        random_seed=seed,
    )
    replica_window_records.append(records)

ensemble = evb.analyze_evb_profile_replicas(
    replica_window_records,
    diabatic_offsets_kcal_mol=(0.0, 0.0),
    coupling_matrix_kcal_mol=((0.0, 5.0), (5.0, 0.0)),
    n_gap_bins=60,
    min_points_per_bin=20,
    minimum_replicas_per_bin=2,
    input_energy_unit="kJ/mol",
    output_directory="results/replicas",
)

print(ensemble["aggregate"]["thermochemistry"]["dg_dagger_mean_kcal_mol"])
print(ensemble["aggregate"]["thermochemistry"]["dg0_mean_kcal_mol"])
```

Replica analysis writes per-replica folders (`rep_000`, `rep_001`, ...) and ensemble files:
- `evb_replica_profile_ensemble.csv`
- `evb_replica_thermochemistry.txt`

### 3. Replace Reactive LJ Walls with a Soft Exponential

For bond-forming/bond-breaking atom pairs (detected from partial bond terms across EVB states),
you can derive soft-exponential parameters and add a correction force that replaces LJ only for
those reactive pairs:

```python
from evbopenmm import EVBSystem

evb = EVBSystem(reactant="state_a.pdb", product="state_b.pdb")

# Optional: inspect which canonical atom pairs are classified as forming/breaking.
reactive_pairs = evb.get_forming_breaking_bond_pairs()

# Derive A, beta for U_exp(r) = A*exp(-beta*r), matched to LJ energy/slope at r_match.
params = evb.derive_reactive_soft_exponential_parameters(
    match_distance_scale=0.85,
    beta_bounds_per_nm=(5.0, 80.0),
)

# Add lambda-scaled correction:
# lambda * (A*exp(-beta*r) - LJ), so LJ is replaced by the soft wall for reactive pairs.
soft_forces = evb.add_reactive_soft_exponential_vdw_forces()
```

Parameter derivation workflow used by `derive_reactive_soft_exponential_parameters(...)`:
- Identify only bonds present in a subset of states (`forming/breaking`).
- Read each selected pair's state-specific LJ `(sigma, epsilon)` from the EVB VdW bond terms.
- Choose `r_match = match_distance_scale * sigma` (kept in the repulsive LJ region).
- Start from the LJ-matched `(A, beta)` at `r_match`.
- If the pair is bonded in another EVB state, use that bond `r0` as a guard distance and enforce:
  - reduced soft slope at guard distance (`slope_ratio_max`, default `0.2`),
  - optional soft-energy cap at guard distance (`guard_energy_cap_kcal_mol`, default `2.0`).
- If needed, downscale the match energy and re-solve `beta` so guard constraints are satisfied.

### 4. Plot LJ vs Soft-Exponential for One Atom Pair

```python
pair = (("reactant", "A", 1, "C1"), ("reactant", "A", 1, "CL1"))

plot_data = evb.plot_pair_lj_soft_exponential_curves(
    pair,
    states=("reactant", "product"),
    energy_unit="kcal/mol",
    n_points=500,
    derive_if_missing=True,
    title="Reactive C–Cl pair",
)

# Optional:
# plot_data["figure"].savefig("pair_lj_soft.png", dpi=200)
```

### 5. Define QM Region for Mixed Nonbonded Treatment

```python
from evbopenmm import EVBSystem

evb = EVBSystem(reactant="state_a.pdb", product="state_b.pdb", auto_build_forces=False)

# Atom refs can come from any state; they are canonicalized internally.
evb.set_qm_atoms(
    [
        ("reactant", "A", 12, "C1"),
        ("reactant", "A", 12, "O1"),
        ("product", "A", 15, "H1"),
    ]
)

# Build nonbonded forces with:
# - MM-MM pairs at cutoff
# - QM-QM and QM-MM pairs with no cutoff
evb.add_vdw_forces(mm_cutoff_distance=1.2)  # nm
```

Available helpers:
- `set_qm_atoms(atom_refs)`
- `get_qm_atoms()`
- `clear_qm_atoms()`

### 6. Solvate All EVB States Consistently (for SCAAS)

```python
from evbopenmm import EVBSystem

evb = EVBSystem(reactant="state_a.pdb", product="state_b.pdb")

solvation = evb.solvate_states(
    reference_state="reactant",
    model="tip3p",
    padding=1.0,          # nm
    ionic_strength=0.15,  # mol/L
    neutralize=True,
)

print(solvation["added_atom_count_per_state"])
print(solvation["state_atom_counts"])
```

`solvate_states(...)` builds solvent once on the reference state and appends the same solvent block to every EVB state, preserving atom ordering needed by EVB mapping and SCAAS water indexing.

### 6b. Solvate as a SCAAS-Matched Sphere (Safe v1)

```python
from evbopenmm import EVBSystem

evb = EVBSystem(reactant="state_a.pdb", product="state_b.pdb")

solvation = evb.solvate_scaas_sphere_states(
    solvent_radius_angstrom=18.0,
    sphere_buffer_angstrom=2.0,
    neutralize=False,      # required in v1
    ionic_strength=None,   # required in v1
)

print(solvation["solvation_mode"])
print(solvation["sphere_center_nm"])
print(solvation["added_atom_count_per_state"])
```

`solvate_scaas_sphere_states(...)` builds a temporary cubic solvent box, trims it to the requested sphere radius, and then clears periodic box vectors on the final EVB topologies.

### 7. Enable Q6-Style SCAAS Sphere Restraints

```python
from evbopenmm import EVBSystem

evb = EVBSystem(reactant="state_a.pdb", product="state_b.pdb")

scaas_info = evb.enable_scaas(
    solvent_radius_angstrom=18.0,
    center=(0.0, 0.0, 0.0),  # nm; omit to auto-center on non-solvent atoms
    polarisation=True,
    charge_correction=True,
    pol_update_steps=100,  # Q6 default (for 1 fs step this is 100 fs)
    backend="cpp",
)

print(scaas_info["solvent_molecule_count"])
print(scaas_info["shells"])
```

Current backend policy:
- The supported SCAAS path is the native `backend=\"cpp\"` implementation.
- `backend=\"cpp\"` is the only accepted public backend value; `python` and `auto` have been removed from the public API.
- OpenMM custom forces provide the radial solvent restraint and shell/excluded positional restraints.
- OpenMM custom compound force provides the polarization term with native shell sorting/target updates on the `cpp` path.
- On CUDA, the conservative host-staged native path remains the default; the faster hybrid device-update path stays
  experimental behind `EVBOPENMM_SCAAS_EXPERIMENTAL_DEVICE_UPDATE=1`.

#### SCAAS native plugin — environment requirement

The `backend="cpp"` plugin (`plugins/scaas/`) has to be **built in the same
environment that Python imports `openmm` from**. The plugin links against
`libOpenMM.so` and (optionally) `libOpenMMCUDA.so`; if those share the C++11
`std::string` ABI with the `openmm` Python module, the plugin loads cleanly.

Mismatches typically surface as a load-time
`undefined symbol: _ZNK6OpenMM18CustomTorsionForce17getEnergyFunctionB5cxx11Ev`
error — that symbol tag (`B5cxx11`) indicates the new C++11 `std::string` ABI,
and the mismatched library carries the old one (no suffix). This happens when
the plugin is built under a conda `openmm` env (new ABI) but Python imports
the `OpenMM.libs` pip-wheel from the base env (old ABI), or vice versa.

The supported workflow is to build and run SCAAS from a single conda env that
has:

- OpenMM (headers + libraries)
- (optional) `libOpenMMCUDA.so` if you intend to run the CUDA platform
- Python with `openmm` importing from that env's site-packages

Check backend availability programmatically with
`evbopenmm.scaas.is_scaas_cpp_backend_available()`. The SCAAS test suite is
decorated with `tests/_backend_markers.py::requires_scaas_cpp_backend` and
skips gracefully when the plugin is unavailable.

### Notes

- The EVB profile APIs are N-state-capable; the examples above show the two-state case. For N ≥ 3 states, build an N-entry `diabatic_offsets_kcal_mol` and an N×N `coupling_matrix_kcal_mol` (symmetric, zero diagonal) and pass `reaction_coordinate_states=(i, j)` to name the two states defining the gap.
- `collect_mapping_fep_samples(...)` stores energies in `kJ/mol`; pass `input_energy_unit="kJ/mol"` when analyzing those records directly.
- For direct external data analysis, use `analyze_evb_profile_records(...)` and `aggregate_evb_profiles(...)` from `evbopenmm.analysis`.

---

## References

1. Warshel, A. & Weiss, R. M. *J. Am. Chem. Soc.* **102**, 6218–6226 (1980).  
2. Åqvist, J. & Warshel, A. *Chem. Rev.* **93**, 2523–2544 (1993).  
3. Kamerlin, S. C. L. & Warshel, A. *WIREs Comput. Mol. Sci.* **1**, 30–45 (2011).

---

© 2025 EVBOpenMM Development Team – GNU License
