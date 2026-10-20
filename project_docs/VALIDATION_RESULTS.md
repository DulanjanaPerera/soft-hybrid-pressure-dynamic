# Validation and benchmark results

These are previously recorded results, not rerun during migration. Full configurations, tables, source links, and limitations remain in the dated [history](PROJECT_LOG.md). MATLAB version was R2025a unless stated otherwise.

| Date / check | Configuration and expected behavior | Recorded result | Interpretation / limit |
| --- | --- | --- | --- |
| Oct 5: local interfaces | validate_pressure_local_interfaces; p=[0,2000,1000] Pa, L=[0.01;0.278;0.015] m, r=0.013 m, K=0.0226, A=pi*(0.013/2)^2 m²; finite dimensions and straight limits | Passed; area sensitivity and basic geometry checked | Callable interfaces only |
| Oct 5: derivatives/integrals | validate_pressure_local_models; four poses through 10 kPa; centered steps 2 Pa / 20 Pa; 24-node quadrature | Passed; largest derivative scaled error 4.529e-6; largest quadrature error 2.914e-7; S derivative identity error 4.614e-16 | Selected-pose numerical agreement; not an operating envelope |
| Oct 5: recursive core | validate_armS_pressure_core; lengths [0.278,0.260,0.290] m, masses [0.10,0.08,0.12] kg; three poses, 20-node direct integration | Passed; dM error <=1.509e-6; mass quadrature error <=2.298e-8; gravity absolute difference <=1.234e-12; all Cholesky checks passed | Translational distributed-mass model only; max rotation orthogonality error 8.803e-8 |
| Oct 5: passive dynamics | validate_armS_pressure_dynamics; 2 s, 401 samples, Q=0, D=8e-11 I, RelTol=1e-8 | Passed; E: 0.928136573 to 0.745525065 J; energy-balance residual 4.758e-5 of energy drop; elastic legacy parity <=9.115e-10; potential-gradient error <=7.337e-8 | Internal consistency; max coordinate 13642.2 Pa exceeds prior 10 kPa local check; no pressure bounds |
| Oct 5: geometry/video | validate_armS_pressure_geometry; straight and trajectory poses; altered sensor offsets; base yaw 25°, pitch 35°, translation | Passed; joint gap and sensor-offset effect zero; world covariance error 4.025e-16; 51-frame video decoded | Mathematical geometry/rendering, not sensor mounting calibration |
| Oct 6: length/pressure geometry | Three poses, four xi values; ell2,ell3=(3Ar²/(2K))*(p2,p3) | Max position/rotation norm discrepancy 9.458e-10; matching example angle about 0.547622952 rad | Kinematic parity only; full dynamics/sign equivalence unverified |
| Oct 7: MEX RHS | build_armS_pressure_mex; four checked configurations | Accepted at 1e-4 tolerance; max scaled RHS component error 3.140e-5 | Approximate counterpart; no bound outside checked poses |
| Oct 7: MEX trajectory | 2 s passive ode15s; 401 samples, RelTol=1e-8 | Max coordinate difference 0.0271519 Pa, rate difference 0.568538 Pa/s | Symbolic/MEX parity on one trajectory |
| Oct 7: performance | Initial-state RHS benchmark; separate 10 s / 2001-sample runner | RHS 0.00301123 s symbolic vs 0.000192656 s MEX; 10 s run integrated in 0.985 s, median MEX RHS 0.000188682 s | Host-specific; excludes rendering; throughput is not closed-loop control frequency |
| Oct 13: angle inverse | validate_config2pressure_reference_N3; unequal parameters, straight pose, central difference step 1e-4 s | Passed; plane-angle error 2.220e-16 rad; rate discrepancy 2.141e-6 Pa/s; acceleration discrepancy 8.822e-5 Pa/s² | Algebra/code checks, not pneumatic calibration |
| Oct 19: trial export helper | save_robosoft_trial smoke test | Runtime verification blocked: initial preferences failure, then license error 63 (-63.2) | Must verify in licensed MATLAB before experimental use |

## Experimental evidence
No hardware identification, held-out command-to-motion validation, or controlled tracking results are recorded here. The manuscript compilation recorded on October 19 produced a four-page preliminary PDF with resolved references and one nonblocking underfull bibliography warning; this is document verification, not model evidence.

## 2026-10-20 documentation verification
Migration checks preserve all original historical text apart from relative-link relocation, confirm all seven documentation files exist, and check local Markdown destinations. No MATLAB or Simulink execution is claimed.
