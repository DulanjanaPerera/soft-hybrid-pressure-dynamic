# Current model and software architecture

## Immediate experimental system
One section, three PMAs, two modeled bending degrees of freedom.

Commanded regulator setpoints p_cmd (3) -> unknown regulator/PMA response -> effective mechanical effort Q (2) -> mechanical state [q_eq;dq_eq] (4) -> theta,phi and tip XYZ.

Actual chamber pressure p_ch is unmeasured. Existing NDI.slx and Pressure_actuationOnly.slx and support scripts are acquisition candidates; their interfaces have not yet been established. Do not equate simulation timestamps with sensor acquisition timestamps without verification.

Single-section baseline: [runner](../runDynamicSimulation_pressure_2DoF_MSF.m), [RHS](../arm_dynamics_pressure_2DoF_MSF.m), M_pressure_2DoF_MSF.m, C_pressure_2DoF_MSF.m, G_and_K_pressure_2DoF_MSF.m. Legacy comments require checking against the coordinate/input distinctions below.

## Coordinate and input contracts
- q_eq: two Pa-valued pressure-equivalent shape coordinates per section, not pressure sensor readings.
- Three-section ordering: [p12;p13;p22;p23;p32;p33]; local generated input [0,p2,p3] is a row with a reference slot.
- Q (legacy single-section tau): generalized effort in J/Pa.
- p_cmd: commanded regulator pressure setpoints; p_ch: distinct actual gauge pressures.
- theta in pressure2config is bending-plane direction; phi is bend magnitude, with legacy wrapping. Confirm experimental conventions.
- L(:,i)=[base sensor offset; flexible length; tip sensor offset]. Physical backbone and mass use [0;L(2,i);0].
- g is base-frame acceleration; world rendering uses Rbase and base position.

## Three-section mechanical implementation
M(q) qdd + C(q,dq)dq + G(q) + F_el(q) + D dq = Q.

| Responsibility | Modules |
| --- | --- |
| Local transform and pressure derivatives | HTM_nume_mex.m, LocalJacob_nume_mex.m |
| Distributed position/derivative moments | integratedPosition_nume.m, integratedPositionProduct_nume.m, integratedPositionDerivativeProduct_nume.m, integratedJacobianProduct_nume.m |
| Symbolic/interpreted recursion | armS_pressure_core.m: M,C,G,dM; armS_pressure_dynamics.m: state RHS |
| Elasticity and energy | pressure_elastic_force.m, armS_pressure_energy.m |
| Compiled path | armS_pressure_entry.m, armS_pressure_core_mex.m, pressure_local_taylor.m, pressure_local_integrals_quadrature.m, build_armS_pressure_mex.m, armS_pressure_mex.mexw64 |
| Run and view | runDynamicSimulation_pressure_N3.m, simulate_armS_pressure_passive.m, armS_pressure_geometry.m, animate_armS_pressure.m, plot_armS_pressure_history.m |
| Reference generation | config2pressure_reference_N3.m; legacy pressure2config/pressure2length helpers |
| Checks | validate_pressure_local_*.m, validate_armS_pressure_*.m, validate_config2pressure_reference_N3.m |

Mass is uniform in normalized xi; only translational kinetic energy is included. C follows the Christoffel construction from dM. G is the gradient of V=-sum(m integral(g' x) dxi). D is symmetric positive semidefinite Rayleigh damping. No pneumatic states or validated physical command allocation are implemented in this three-section core.

## Experiment data flow
Retain raw command and pose streams with original timestamps and metadata -> verify clock/units/frame conventions -> documented alignment/filtering in analysis copies -> trial-level train/validation split -> bounded forward-simulation fitting -> held-out angle and XYZ comparisons -> conditional feedback tracking comparison.

See [logging contract](../ROBOSOFT_2027_DATA_LOG_SPEC.md) and [export helper](../save_robosoft_trial.m). MATLAB is required for runtime checks; MEX build additionally uses MATLAB Coder and MSVC. Export helper runtime testing is unresolved.

## Detailed references
The [Cosserat handoff](../cosserat_rod_pressure_model_handoff.md) preserves detailed pressure/length contracts and proposed rod development. It does not establish an implemented Cosserat model. The [original recursion handoff](../recursive_dynamics_pressure_coordinate_handoff.md) is historical context. Current state and newer decisions take precedence over superseded plans.
