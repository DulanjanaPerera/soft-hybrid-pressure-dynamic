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

## 2026-10-05 — Step 3: interpreted three-section recursive core

**Implementation:** Added [armS_pressure_core.m](armS_pressure_core.m), returning `M` (6×6), `C` (6×6), `G` (6×1), and `dM` (6×6×6) for `q=[p12;p13;p22;p23;p32;p33]` (Pa) and `dq` (Pa/s). Each column of `params.L` is one section's `[base sensor offset; flexible section length; tip sensor offset]` in meters; the validation uses three unequal section lengths. Section mass is uniform in normalized material coordinate, `dm_i=m_i dxi`. Only translational kinetic energy is modeled. Each section's mass moments and the tip transform used for propagation receive `Lbody=[0;params.L(2,i);0]`, so sensor offsets affect neither the mass terms nor the physical connection between section backbones. The local pressure exports are used directly; no length-coordinate transformation is applied.

**Dynamics convention:** `params.g` is gravitational acceleration in the model frame. The potential is `V(q)=-Σ_i m_i∫gᵀx_i(q,xi) dxi`, and `G=∂V/∂q`. The intended equation is `M(q) qdd + C(q,dq) dq + G(q) = Q`, with the pressure-coordinate generalized input `Q` still to be derived. `C` is assembled from `dM` with the Christoffel convention used in the handoff. The current local-local mass block treats the polynomial `R` as orthogonal; the approximation is measured below.

**Independent validation:** Added [validate_armS_pressure_core.m](validate_armS_pressure_core.m). The reference composes section transforms directly, estimates global point Jacobians by centered 2 Pa pressure differences, and integrates `JᵀJ` and `-Jᵀg` with 20-node Gauss–Legendre quadrature. It checks `dM` with centered 20 Pa differences of `M`, `Mdot-2C` skew symmetry, mass symmetry and Cholesky positive definiteness, and invariance when only the massless sensor offsets change. At the second pose, it also compares `G` with a finite-difference gradient of the directly integrated potential.

The test parameters were `L(2,:)=[0.278,0.260,0.290]` m, `r=[0.013;0.0125;0.014]` m, `K=[0.0226;0.021;0.024]` (current geometry parameters; units remain to be reviewed), `A=pi*(r/2).^2` m², `m=[0.10;0.08;0.12]` kg, `g=[0;0;-9.81]` m/s², and `dq=[100;-70;50;20;-40;80]` Pa/s. Sensor offsets were nonzero and distinct in the input `L` matrix.

**MATLAB R2025a results:** `matlab -batch "validate_armS_pressure_core"` passed.

| Pressure pose `q` (Pa) | M symmetry | `dM` vs finite difference | `Mdot-2C` skew | M vs direct quadrature | G absolute difference | Max `RᵀR-I` Frobenius norm |
|---|---:|---:|---:|---:|---:|---:|
| `[0;0;0;0;0;0]` | `0` | `0` | `0` | `1.928e-8` | `0` | `0` |
| `[4500;1500;3000;1000;2000;4000]` | `2.991e-17` | `1.509e-6` | `9.763e-17` | `1.998e-8` | `1.136e-12` | `3.388e-10` |
| `[10000;0;2000;6000;5000;1000]` | `0` | `1.401e-6` | `1.041e-16` | `2.298e-8` | `1.234e-12` | `8.803e-8` |

All three mass matrices passed Cholesky. The gravity-to-potential-gradient scaled error at the second pose was `1.020e-6`. Sensor-offset-only changes left `M`, `C`, `G`, and `dM` identical in all three poses. Ratios in the table use the norms and floors in the validation script; the G column is an absolute Euclidean difference.

**Limits:** This validates the interpreted distributed-mass core at three chosen poses. It does not establish accuracy at larger bends, physical actuation and elastic/damping laws, time integration, animation, or MEX parity. The rotational kinetic energy and any sensor hardware mass remain excluded by model choice; the polynomial rotation's nonorthogonality should be monitored when expanding the operating range.

## 2026-10-05 — Step 4: passive pressure-coordinate dynamics

**Elastic law and provenance:** Added [pressure_elastic_force.m](pressure_elastic_force.m). The elastic part of the existing single-section [G_and_K_pressure_2DoF_MSF.m](G_and_K_pressure_2DoF_MSF.m) at `g=0` is linear in the two independent pressures. For section `i`, let `a_i=A_i²r_i²/K_i` and `b_i=9r_i²/(4K_i)`. The implemented block is `H_i=a_i*[3+b_i*k(2,i), 3/2; 3/2, 3+b_i*k(3,i)]`, with `U_i=1/2 q_iᵀH_i q_i` and `F_elastic,i=H_i q_i`. `k(1,i)` does not appear explicitly in this block; the example selects `K_i=(3/20)k(1,i)L(2,i)r_i²`, matching the legacy one-section parameter choice. This is a pressure-coordinate constitutive model inherited from the existing expression, not a new pneumatic work/actuation derivation. Its units and physical calibration still need review.

**Equation and implementation:** Added [armS_pressure_dynamics.m](armS_pressure_dynamics.m) for state `X=[q;dq]` and `M qdd+C dq+G+F_elastic+D dq=Q`. `D` is a symmetric positive semidefinite `6×6` Rayleigh damping matrix; `Q` is an explicitly supplied work-conjugate generalized force in J/Pa. No pressure-command, valve, or regulator dynamics have been inferred. Added [armS_pressure_energy.m](armS_pressure_energy.m), calculating `E=1/2 dqᵀM dq+ΣU_i-Σm_i gᵀ(P_base,i+R_base,i μ_i)`. The gravitational potential uses each section's offset-free backbone first moment and the physical tip transform. The sensor offsets remain massless.

**Interpreted example:** Added [simulate_armS_pressure_passive.m](simulate_armS_pressure_passive.m). It uses `L(2,:)=[0.278,0.260,0.290]` m, `r=[0.013;0.0125;0.014]` m, `k=3200 N/m` for each PMA, `A_i=π(r_i/2)²` m², `m=[0.10;0.08;0.12]` kg, `g=[0;0;-9.81]` m/s², `D=8e-11 I`, `Q=0`, `q(0)=[10000;0;2000;1000;1000;0]` Pa, and `dq(0)=0`. `ode15s` integrates for 2 s with 401 reported samples, relative tolerance `1e-8`, and absolute tolerances `1e-5` Pa for `q` and `1e-4` Pa/s for `dq`.

**MATLAB R2025a validation:** [validate_armS_pressure_dynamics.m](validate_armS_pressure_dynamics.m) passed with `matlab -batch "validate_armS_pressure_dynamics"`. At initial, midpoint, and final poses, the largest relative difference between the new elastic block and the legacy `g=0` function was `9.115e-10`; the largest relative error between a centered 2 Pa finite difference of `U+V` and `F_elastic+G` was `7.337e-8`. The elastic matrix was symmetric positive definite. The passive run remained finite, energy fell at every sampled interval, and `E(0)=0.928136573 J`, `E(2)=0.745525065 J`. The residual `E(2)-E(0)+∫dqᵀD dq dt`, using trapezoidal integration over reported samples, was `4.758e-5` of the energy drop. Maximum absolute pressure coordinate was `13642.2` Pa; the largest sampled energy increment was `-2.335e-5` J.

**Limits:** The largest trajectory coordinate exceeds the 10 kPa pose used in the previous local-model checks; this simulation demonstrates numerical integration and energy consistency, not a calibrated operating envelope. Pressure coordinates can change sign during passive release; no pressure bounds or actuator constraints are imposed. The elastic law, `K_i` choice, and damping value are inherited or illustrative and require experimental calibration. Animation, explicit actuation input mapping, rotational inertia, and MEX generation remain future steps.

## 2026-10-05 — Step 5: three-section trajectory animation

**Geometry:** Added [armS_pressure_geometry.m](armS_pressure_geometry.m). For every reported material coordinate `xi`, it composes each section's pressure-local position and rotation with the preceding physical tip transforms. The renderer uses `Lbody=[0;L(2,n);0]` for all three backbones and propagates at `xi=1`. The two sensor offsets in each `params.L(:,n)` do not add tube length or separate adjacent sections. The section frame's first two rotation columns orient each tube cross-section. This visualization does not claim a sensor measurement frame or hardware geometry; those require separate mounting information.

**Animation:** Added [animate_armS_pressure.m](animate_armS_pressure.m). It draws the three sections as distinct colored tubes, marks the tip, and uses fixed axes fitted to the entire selected trajectory. Options include frame stride, mesh resolution, tube radius, and video output. The passive simulation from Step 4 was rendered with `frameStride=8`, 31 backbone samples per section, 16 circumferential samples, radius `0.015 m`, and 25 frames/s. Saved [results/armS_pressure_passive.mp4](results/armS_pressure_passive.mp4) and a first-frame [preview](results/armS_pressure_preview.png). The MP4 is 2.040 s at 1126×878 pixels; the final simulation sample is included even when the stride does not land on it.

**MATLAB R2025a checks:** [validate_armS_pressure_geometry.m](validate_armS_pressure_geometry.m) tested straight, initial, midpoint, and final trajectory poses. The largest section joint gap was `0 m`; changing all sensor offsets by `0.1 m` changed no rendered position, orientation, or base point (`0` maximum difference). A three-frame headless rendering passed. `VideoReader` decoded all 51 frames of the produced MP4. The first and final frames were inspected visually for connected colored sections and correct camera framing.

**Limits:** This is a visualization of the interpreted mathematical trajectory. It does not establish the physical pressure range, actuator input law, sensor marker placement, wall deformation, or correspondence with an experiment. MEX generation remains the next implementation step.

## 2026-10-05 — Editable three-section simulation and video script

**Request and interface:** Added [runDynamicSimulation_pressure_N3.m](runDynamicSimulation_pressure_N3.m), following the user-editable layout of the repository's one-section `runDynamicSimulation_pressure_2DoF_MSF.m`. The script exposes each section's `L(:,i)`, `r_i`, `A_i`, `m_i`, PMA stiffness column `k(:,i)`, bending stiffness `K_i`, gravity, the 6×6 pressure-coordinate damping matrix `D`, six initial pressures, six initial rates, simulation time and solver tolerances, and the animation/video settings. The default `K_i=(3/20)k(1,i)L(2,i)r_i²` is the inherited example relation and can be replaced with a measured 3×1 vector. The script uses the interpreted [armS_pressure_dynamics.m](armS_pressure_dynamics.m) and [animate_armS_pressure.m](animate_armS_pressure.m) without altering their equations.

**Input convention:** `inputForce` may be a constant 6×1 vector or a function handle `@(time,state)` returning a 6×1 vector. Its entries are generalized forces conjugate to the six pressure coordinates, in J/Pa. It is not a chamber pressure command, valve signal, or regulator input. The default is zero, giving the Step 4 passive release. `t`, `X`, `params`, and `videoFile` remain in the MATLAB workspace after the script finishes.

**Run check:** In MATLAB R2025a, `matlab -batch "runDynamicSimulation_pressure_N3"` completed the default simulation with 401 samples over 2.000 s (17.60 s elapsed on this run). It recorded [results/armS_pressure_run.mp4](results/armS_pressure_run.mp4); `VideoReader` reopened it and reported a 2.040 s video. The output is ignored by Git because rerunning the script with changed parameters replaces it. The script itself and this log are versioned.

**Limit:** This is a configurable numerical experiment using the present pressure-coordinate elastic/damping assumptions. It does not provide a physical mapping from chamber pressure commands to generalized force. MEX generation and calibrated actuation remain future work.

## 2026-10-05 — Rotatable base and reusable figures

**Research decision and implementation:** The three-section runner now exposes base roll, pitch, and yaw angles in degrees and a world-frame base position. It forms `Rbase=Rz(yaw)Ry(pitch)Rx(roll)`. [armS_pressure_geometry.m](armS_pressure_geometry.m) begins recursive rendering from that pose, so every backbone point, section frame, and tip is expressed in world coordinates. The runner stores `params.gWorld` and supplies `params.g=Rbaseᵀ params.gWorld` to the dynamics core, whose generalized coordinates and gravity equations are in the base frame. This makes a physical base tilt alter gravity relative to the arm, while a pure change of base position changes only the constant reference of gravitational potential. Zero angles reproduce the previous orientation. The researcher's existing edits to section lengths, mass, gravity direction, damping, initial state, and ten-second run were preserved.

**Figures:** [animate_armS_pressure.m](animate_armS_pressure.m) now clears and reuses MATLAB Figure 1 on each invocation, rather than creating a new numbered figure. A base marker and trajectory-fitted fixed axes also accommodate a moved base. Added [plot_armS_pressure_history.m](plot_armS_pressure_history.m), which clears and reuses Figure 2 and shows two independent pressure coordinates (kPa) for each of the three modules in separate stacked panels. The runner enables this plot by default with `makePressurePlots=true` and still records the animation when `videoOptions.saveVideo=true`.

**MATLAB R2025a checks:** Extended [validate_armS_pressure_geometry.m](validate_armS_pressure_geometry.m) with a 25° yaw, 35° pitch, and translated base. Across straight, initial, midpoint, and final trajectory poses, the maximum world-transform covariance error was `4.025e-16`; joint gaps and sensor-offset effects remained zero. The researcher's current `runDynamicSimulation_pressure_N3.m` settings completed 2001 samples over 10.000 s, recorded a 10.040 s MP4, and produced Figures 1 and 2 with six pressure traces. Repeated two-frame rendering with a 180° base roll reused the same Figure 1 and Figure 2 handles; Figure 2 still contained six traces.

**Limit:** `params.gWorld` in the runner is converted to base-frame `params.g` before simulation. Direct calls to the core must continue supplying gravity in the base frame. Rotating the base changes orientation and gravity loading, not the constitutive or generalized-input law.

## 2026-10-05 — Interpretation of pressure coordinates and physical actuation

**Research discussion:** A passively falling, vented arm can change shape while its chamber gauge pressures stay approximately zero. Therefore the six simulated coordinates currently called pressures must not be interpreted automatically as measured chamber pressures. They are pressure-valued parameters of the chosen shape map (an equivalent-pressure shape chart). A nonzero value during the passive simulation labels a deformed shape; it does not claim the vented chamber has acquired that pressure. A sealed chamber is different: deformation can change its volume and pressure even with no valve command. Whether the current pressure-to-geometry map remains valid under gravity and external loading must be checked, especially if it was originally fitted or derived from unloaded quasi-static actuation equilibria.

**Meaning of the current input:** Since `q` is measured in Pa, its generalized force satisfies `δW=Q_qᵀδq` and has units `J/Pa=m³`. This is a virtual-work conjugate quantity, not a physical pressure command. If a mechanical shape vector `s` obeys `s=f(q)` with Jacobian `J_f=∂s/∂q`, then a mechanical generalized force `F_s` maps to `Q_q=J_fᵀF_s`. For chamber gauge pressures `p_ch` and volume map `V(s)`, pneumatic virtual work is `δW=p_chᵀδV`, yielding `Q_q=(∂V/∂q)ᵀp_ch=J_fᵀ(∂V/∂s)ᵀp_ch` under the stated holonomic map. These equations are a proposed physical input derivation; chamber volume functions and the validity of the shape map have not been established for this arm.

**Modeling consequence:** For a physically interpretable passive/actuated comparison, keep mechanical shape (length or curvature, or a clearly labeled pressure-equivalent shape coordinate) distinct from actual chamber-pressure states. Vented passive motion uses approximately zero chamber gauge pressure but may have nonzero shape coordinates. A pneumatic model would add chamber volumes, mass-flow/valve dynamics, pressure states, and a virtual-work coupling. The present `inputForce` remains an abstract generalized force until that derivation and calibration are done. The Figure 2 title and runner comments were updated to label the plotted trajectories as pressure-equivalent shape coordinates; they are not predicted sensor readings.

**Research context:** This separation of mechanical coordinates and pneumatic pressure states, and the pressure-volume work coupling, is consistent with [Stölzle and Della Santina, *Piston-Driven Pneumatically-Actuated Soft Robots: Modeling and Backstepping Control*](https://pure.tudelft.nl/ws/portalfiles/portal/103544037/Piston_Driven_Pneumatically_Actuated_Soft_Robots_Modeling_and_Backstepping_Control.pdf) and [*Nonlinear energy-based control of soft continuum pneumatic manipulators*](https://link.springer.com/article/10.1007/s11071-021-06817-1). Their constitutive and pneumatic assumptions are not adopted automatically here.
