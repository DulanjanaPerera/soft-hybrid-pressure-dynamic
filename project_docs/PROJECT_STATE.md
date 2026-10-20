# Project state
Last Updated: 2026-10-20
Current Milestone: Commission single-module command and NDI pose logging for experimental identification

## Project objective
Identify and experimentally validate command-to-motion dynamics of one hybrid continuum module with three PMAs and two bending degrees of freedom, then evaluate model-assisted tracking if supported by the evidence. This supersedes the immediate three-section paper scope. The three-section recursive implementation remains available.

## Current architecture / model
The mechanical model uses pressure-equivalent shape coordinates q_eq, not measured chamber pressures. Generalized effort Q is work-conjugate, in J/Pa. Three regulator pressure setpoints p_cmd drive an uncalibrated regulator/PMA/mechanical chain; actual chamber pressures p_ch are unmeasured. Available outputs are tip XYZ and theta,phi, potentially derived from the same NDI sensing pipeline.

The immediate baseline is [single-section dynamics](../arm_dynamics_pressure_2DoF_MSF.m) and its [runner](../runDynamicSimulation_pressure_2DoF_MSF.m). A separate six-coordinate, three-section model has interpreted and MEX implementations, passive simulation, geometry, and animation. See [architecture](MODEL_ARCHITECTURE.md).

## Current implementation status
- Three-section numerical verification, passive energy checks, animation, and MEX parity/timing are recorded; these are not physical model validation.
- [Paper plan](../ROBOSOFT_2027_SINGLE_SECTION_PAPER_PLAN.md) and [data-log specification](../ROBOSOFT_2027_DATA_LOG_SPEC.md) define the current experiment.
- [save_robosoft_trial.m](../save_robosoft_trial.m) exports timestamped tables and metadata without overwriting raw trial files. MATLAB runtime verification was blocked by licensing error 63 (-63.2), after an initial preferences-directory failure.
- The researcher reports NDI and pressure-command Simulink models. NDI.slx, Pressure_actuationOnly.slx, and NDI support scripts are present; their acquisition paths and timing have not been inspected.
- The historical log records a preliminary manuscript in the sibling Paper repository. Its main branch serves Overleaf; the latest recorded cleanup commit is a772e2a. This migration did not inspect or modify that repository.

## Validated components
Selected-pose local derivatives/integrals, recursive mass/gravity and Christoffel consistency, passive energy behavior, base/world geometry covariance, video decoding, and limited MEX/symbolic parity are documented in [validation results](VALIDATION_RESULTS.md). None establishes command-to-motion prediction on hardware.

## Current task
Inspect the existing Simulink models and support code, identify command and NDI pose channels, timestamps, sample rates, units, calibration, angle conventions, and export paths. Produce and inspect one commissioning trial before identification chirps.

## Open issues
See [active issues](OPEN_ISSUES.md). Key gaps: command/pose synchronization; unmeasured chamber pressure; effective input coupling and parameter identifiability; K definition/units; high-bend numerical approximation; historical gravity-sign mismatch; unverified export helper.

## Important current assumptions
- Fitted stiffness, damping, gain, and any lag are effective command-to-motion parameters unless independent measurements justify stronger interpretation.
- The three physical commands are feasible nonnegative pressure setpoints; the shape chart's zero-sum relation is not a physical pressure constraint.
- XYZ and angles from the same tracker are correlated outputs. Verify theta/phi definitions before comparing with model helpers.
- Existing numerical tests cover listed configurations only. Massless sensor offsets and exclusion of rotary inertia define the current three-section model.
- No experimental results, physical controller validation, or calibrated regulator dynamics are established in this documentation.

## Next steps
1. Inspect NDI.slx and Pressure_actuationOnly.slx without operating hardware.
2. Verify pose frames, calibration, angle definitions, command scaling, clock origin, rates, missing samples, and logging.
3. Collect a modest commissioning trial with synchronized command and pose records and metadata.
4. Freeze training/validation trial splits, collect informative excitations, and fit forward-simulated angles with bounded effective parameters.
5. Compare held-out predictions against static and low-order dynamic baselines; proceed to control only if justified.

## Documentation maintenance
Start sessions here; check this date against newer relevant work. Read the specialized document needed for the task, and search historical entries only when necessary. After significant work append [history](PROJECT_LOG.md), update this state/date and affected specialized documents. History is append-only after migration. Preserve uncertainties and distinguish numerical checks from experimental evidence.
