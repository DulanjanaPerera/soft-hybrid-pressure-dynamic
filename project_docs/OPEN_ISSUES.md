# Active open issues

## Acquisition and synchronized logging
Evidence: NDI.slx and Pressure_actuationOnly.slx and support code are present; the researcher reports using these for tracking/commands. Status: uninspected. Timing, actual sent signals, frame calibration, units, angle conventions, sensor latency, stale samples, and synchronization are unknown. Next: inspect signal paths and exports, then one commissioning trial using the [data-log specification](../ROBOSOFT_2027_DATA_LOG_SPEC.md). No hardware operation has been attempted in this migration.

## Export helper runtime verification
Attempted MATLAB smoke test failed first on preferences permissions, then on licensing error 63 (-63.2). Status: unverified. Next: run save_robosoft_trial in licensed MATLAB with synthetic tables and check schema, metadata, timestamp rejection, and overwrite protection.

## Effective input coupling and unavailable chamber pressures
Only command setpoints and pose are available; no p_ch sensor channel. No calibrated physical Q=B(q)p_ch relation is established. Proposed volume/virtual-work and constant-area formulas are hypotheses, not identified laws. Status: physical bridge unresolved; immediate study is command-to-motion. Next: select a parsimonious effective coupling and validate held-out trajectories. Do not assign a fitted lag solely to the regulator.

## Parameter identifiability and K units
Geometry angle formula requires K in N·m, while inherited K=(3/20)k L r² has N·m² if k is N/m; old N/rad comments also conflict. No resolution implemented. Static data identify actuator-gain/stiffness ratio, not both independently; chart also depends on K. Next: revisit derivation/calibration, fix independently known quantities, inspect parameter sensitivity, and report effective values where necessary.

## High-bend approximation and compiled assembly maintenance
Polynomial rotations lose orthogonality; symbolic integrals and transform quadrature disagree increasingly at high bend. MEX uses compact approximations and a separate recursive assembly. Attempt: four-pose parity and one trajectory comparison passed their recorded criteria. Status: expanded operating envelope unverified. Next: test experimentally relevant poses, monitor errors, and synchronize both assemblies when changing equations.

## Historical length-model gravity sign
Length reference assembles +sum(m J' g), pressure core uses -sum(m J' g); both RHS functions subtract G. Kinematic parity does not settle dynamic equivalence. Status: flagged, not corrected. Next: independent potential-gradient and physical falling-direction checks before reusing or comparing length dynamics.

## Geometry assumptions and sensor interpretation
Massless offsets and connected backbone geometry are numerically checked, but actual sensor mounting frames, hardware mass, extension/common-pressure effects, payloads, and shape-chart validity under load need evidence. Next: compare NDI frames and measured geometry against model assumptions; test common command separately.

## Experimental evidence and paper readiness
No commissioning dataset, identified parameters, held-out model comparison, or hardware control result is established. Historical paper dates/schedule remain a plan, not verified completion or current venue policy. Next: establish usable command/pose data, then static and dynamic baselines and conditional control. Recheck venue requirements when preparing submission.
