# Restoration validation — 17 November 2026

## Newly executed checks

- **21 Python tests pass.** They cover Q topology parsing, native LJ coefficient conversion, state charge/constraints, minimization, transactional restart, frame-wise reconstruction and the portable CLI.
- **Browser checks pass.** The page exposes full-enzyme commands, contains the specific Q6/OpenMM comparison, has no tested mobile horizontal overflow, and keeps parameter controls separate from the fixed paper reference.
- **Installed bundled package checked.** A no-index local pip installation was used to construct all 6,330 atoms. Pure-state and mixed-state energies/forces were finite, the checked coordinate passed a finite-difference test and mapping-linearity residual was zero. This is not full native-Q6 energy/force parity.
- **New end-to-end full-system smoke completed.** Des27.7-IN, 6,330 atoms, 771 waters, CPU with two threads. Bounded descent/L-BFGS relaxation, heating, endpoint equilibration and three mapping windows completed. It performed 480 MD steps at 0.25 fs (0.12 ps) and saved five samples/window, plus PDB, DCD and energy traces. Diagnostics and profile export ran on those records. This is not adequate sampling or a physical barrier measurement. The minimizer did not reach its force-tolerance target; the constraint-normal forces and bounded minimization limits are reported.
- **Historical full-system result reproduced.** The supplied 21-window raw record matches SHA256 `99946c99cd0d8769044d70a611146f4a5438001e728cbd0e4e05a833401978a1`. With 18 bins/minimum 20/window-bin, the reconstructed conditional barrier is 9.572629012215984 kcal/mol. This is reanalysis of the archived run, not a new production simulation or a target fit. Explicitly supplying the same parameters gives the same frame-wise reconstruction.

The smoke profile-export check used three bins and minimum one sample/window-bin solely to test output plumbing. Those settings are not recommended for interpretation. No smoke-derived barrier is claimed.

## Still unvalidated

Native Q6 LRF/SCAAS, restraint and full Hamiltonian parity; independent calibration of transferred coupling/offset and reactive terms; adequate equilibration and timestep stability; overlapping windows, effective sample sizes, independent replicas and a converged activation/reaction free energy. The short and long presets are not certified scientific protocols.

Small test reports are retained in `docs/validation/`. Generated full Systems, trajectories, restart databases and working output directories are not committed. Restoring source and demonstrating execution do not establish exact paper reproduction.
