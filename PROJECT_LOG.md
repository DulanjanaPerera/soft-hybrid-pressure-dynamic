# Pressure-coordinate recursive dynamics project log

This log records modeling decisions, implementation changes, validation evidence, and unresolved questions for the three-section soft continuum arm. Append a dated entry for every code change or material research discussion. Keep measured results and verified calculations separate from assumptions or proposed methods. Record parameters, units, test conditions, and source files so results can be traced when preparing the paper.

## 2026-10-02 — Repository baseline and work branch

- The pressure-dynamics repository is `DulanjanaPerera/soft-hybrid-pressure-dynamic`.
- The then-pending simulation, Simulink, animation, result video, and generated cache files were committed on `main` as `1cbbf3c` (`Save current pressure dynamics and animation work`) and pushed.
- Branch `codex/recursive-dynamic-modeling` was created from that commit and published. Recursive modeling work belongs on this branch.

## 2026-10-05 — Handoff and model scope

**Sources reviewed:** [recursive_dynamics_pressure_coordinate_handoff.md](recursive_dynamics_pressure_coordinate_handoff.md), the local pressure-space `*_nume.m` exports, the historical length-coordinate project branches (`Ozi` and `compare-pointmass-standard`), and the clean `standard-arm-floating-base` implementation.

**Decisions from the researcher:**

- Use pressure itself as the generalized coordinate. The local symbolic transform, Jacobians, Hessians, and section integrals have already been differentiated with respect to pressure. No length-to-pressure coordinate transformation is required.
- Use two independent pressures per section, as in the two-length-change model. For three sections, the intended generalized coordinate ordering is `q = [p12; p13; p22; p23; p32; p33]` in Pa. A local generated function receives `[0, p2, p3]` as a `1x3` row; the leading zero is a reference slot, not a third independent coordinate.
- In each section's local `L = [L(1); L(2); L(3)]` (m), `L(2)` is the flexible section length. `L(1)` and `L(3)` are sensor offsets. They are not lengths of other sections and do not carry mass.
- Sensor geometry may use the full local `L`. For distributed backbone mass integrals, use `Lbody = [0; L(2); 0]` so sensor offsets do not contribute mass. The placement of sensor frames in the future three-section animation and measurement model still needs explicit validation.

**Implementation direction:** adapt the clean reference recursion for pressure-local geometry, then validate local derivatives and integrals, recursive mass and gravity terms, velocity terms, interpreted simulation, animation, and finally MEX parity. The historical length-coordinate functions are references for algorithm structure and validation methods, not substitutes for the pressure exports.

**Open modeling decisions:** the pressure-coordinate generalized input/work law, elastic and damping laws, section-specific physical parameters, gravity convention, and any sensor measurement frames have not yet been established for the three-section model. No three-section pressure dynamics are claimed as validated.

## 2026-10-05 — Step 1: local function interfaces

**Change:** The six generated local functions used the symbol `A` in their expressions without accepting or defining it. Added explicit positive scalar `A` (effective PMA area, m²) as the last input to [HTM_nume_mex.m](HTM_nume_mex.m), [LocalJacob_nume_mex.m](LocalJacob_nume_mex.m), [integratedPosition_nume.m](integratedPosition_nume.m), [integratedPositionProduct_nume.m](integratedPositionProduct_nume.m), [integratedPositionDerivativeProduct_nume.m](integratedPositionDerivativeProduct_nume.m), and [integratedJacobianProduct_nume.m](integratedJacobianProduct_nume.m). Documented the pressure row, `L` meanings, and units; corrected stale length-coordinate labels in relevant derivative comments. The generated symbolic expressions were not changed.

**Validation:** Added [validate_pressure_local_interfaces.m](validate_pressure_local_interfaces.m). On MATLAB R2025a, `matlab -batch "validate_pressure_local_interfaces"` passed. The check used `p=[0,2000,1000]` Pa, `L=[0.01;0.278;0.015]` m, `Lbody=[0;0.278;0]` m, `r=0.013` m, `K=0.0226` (the current model parameter; its documented units need review), and `A=pi*(0.013/2)^2` m². It checked output dimensions, finite values, transform output consistency, straight-backbone endpoints and first moment with zero offsets, and sensitivity to `A`.

**Limit:** This establishes callable interfaces and basic geometry only. It does not yet verify analytic pressure derivatives, integral identities, mass properties, or MEX code generation. Those checks are the next step.

## 2026-10-05 — Step 2: local pressure derivatives and section integrals

**Method:** Added [validate_pressure_local_models.m](validate_pressure_local_models.m). At each pose, centered pressure perturbations check the position and rotation Jacobians against [HTM_nume_mex.m](HTM_nume_mex.m) and the Hessians against differences of those Jacobians. A 24-node Gauss–Legendre rule integrates local point position, Jacobian, and Hessian over `xi ∈ [0,1]` and compares the results with the exported `mu`, `S`, `F`, `E` moments and their pressure derivatives. The validation also checks `S_,a = F_a + F_aᵀ`. The numerical reference uses `Lbody=[0;L(2);0]`, since the sensor offsets have no mass; derivative checks use the full `L` to exercise sensor geometry. Centered-difference steps were 2 Pa for first derivatives and 20 Pa for second derivatives.

**Defect and correction:** The expanded transform returned `NaN` in two rotation entries at exactly `p2=p3=0`. Values approaching zero pressure converged to the identity rotation. [HTM_nume_mex.m](HTM_nume_mex.m) now returns the straight-arm limit `R=I`, `P=[0;0;L(1)+xi*L(2)+L(3)]` at that exact configuration, before evaluating the removable `0/0` terms. The symbolic expression used at nonzero pressure was not changed.

**MATLAB R2025a results:** `matlab -batch "validate_pressure_local_models"` passed. The table reports the largest scaled error among the four derivative comparisons and among the nine quadrature comparisons at each pose. The scale is the numerical reference norm with a `1e-12` floor.

| Local pressure `[0,p2,p3]` (Pa) | `L(2)` (m) | `xi` for derivative check | Largest derivative error | Largest quadrature error |
|---|---:|---:|---:|---:|
| `[0,0,0]` | 0.278 | 0.35 | `1.428e-7` | `1.643e-15` |
| `[0,4500,1500]` | 0.260 | 0.67 | `5.237e-7` | `8.745e-10` |
| `[0,2000,8000]` | 0.290 | 1.00 | `4.529e-6` | `1.403e-7` |
| `[0,10000,0]` | 0.278 | 0.50 | `2.920e-7` | `2.914e-7` |

The largest `S_,a = F_a + F_aᵀ` scaled error was `4.614e-16`. The quadrature comparison includes mixed terms in `F_q` and `E_q`; it does not substitute products of separately integrated quantities. The 10 kPa pose matches the order of the current single-section simulation's initial pressure.

**Limit:** These tests establish numerical agreement at the four listed poses and supplied section lengths, not a validated operating envelope. The expanded expressions show increasing numerical disagreement at larger bends. Three-section recursive assembly, independent mass and gravity checks, physical force laws, trajectory validation, and MEX generation remain untested.
