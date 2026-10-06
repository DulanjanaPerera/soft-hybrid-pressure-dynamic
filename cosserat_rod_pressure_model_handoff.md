# Three-section arm: kinematics and dynamics handoff for a Cosserat rod model

- **Prepared:** 6 October 2026
- **Repository:** `DulanjanaPerera/soft-hybrid-pressure-dynamic`, branch `codex/recursive-dynamic-modeling`
- **Source baseline inspected:** commit `722707d` and the working MATLAB files in this directory.
**Purpose:** Give another implementation chat a traceable description of the present arm model and a staged route to an efficient Cosserat rod dynamic model. This is a technical handoff, not a claim that the current pressure-coordinate simulation is a calibrated pneumatic plant.

## 1. Model at a glance

- There are **three serial flexible sections**. Each has **two independent shape coordinates**. The 6-vector used by the present model is
  `q = [p12; p13; p22; p23; p32; p33]` in Pa, with `X=[q; qdot]` (12×1).
- A section-local generated function receives `p_i=[0,p_i2,p_i3]` as a 1×3 row. The leading zero is a reference slot; it is **not** a third independent generalized coordinate or a claim that a physical first chamber always has zero gauge pressure.
- The six `q` entries are best called **pressure-equivalent shape coordinates**. They parameterize geometry. During passive bending they need not equal measured chamber pressures. Actual chamber pressures and valve states are not modeled.
- `params.L(:,i)=[base sensor offset; flexible section length; tip sensor offset]` (m). `L(1,i)` and `L(3,i)` are massless sensor offsets, **not** other sections' lengths. For mass-bearing backbone geometry and propagation, use `Lbody_i=[0;L(2,i);0]`.
- The existing three-section dynamics include distributed **translational** mass, gravity, a pressure-coordinate elastic potential, and viscous damping. Cross-sectional rotary inertia, shear/extension, torsion, sensor hardware mass, payloads, contacts, and pneumatic states are absent.
- The local transforms are expanded polynomial expressions derived in Maple (the historical handoff identifies a 13th-order Taylor transform). They are verified numerically over selected poses, but their rotation matrices are only approximately orthogonal at large bends. The new rod should use an exact `SO(3)`/`SE(3)` representation where possible.

## 2. Parameters, frames, and coordinate meaning

| Quantity | Present code | Meaning |
|---|---|---|
| `q_i=[p_i2;p_i3]` | 2×1 per section | Shape chart with Pa-valued coordinates; not a chamber-pressure state |
| `L_i=L(2,i)` | m | Flexible backbone arc length in the reference configuration |
| `L(1,i),L(3,i)` | m | Sensor offsets; omitted from physical backbone and mass integrals |
| `r_i` | m | Actuator radial offset from backbone |
| `A_i` | m² | Effective actuator area in the pressure-to-geometry expressions |
| `K_i` | code comment says N/rad | Parameter in the generated geometry; **units need resolution**, section 7 |
| `m_i`/`params.mi(i)` | kg | Total distributed mass of section `i` |
| `k(:,i)` | code comment says N/m | Three local PMA stiffness values used by the inherited elastic law |
| `g`/`params.g` | m/s² | Gravity expressed in the **base frame** for the dynamics core |
| `D` | 6×6 | Symmetric positive semidefinite viscous matrix in the `q` chart |
| `Q`/`inputForce` | J/Pa = m³ | Generalized effort conjugate to `q`; **not** a chamber-pressure command |

The user-editable [runDynamicSimulation_pressure_N3.m](runDynamicSimulation_pressure_N3.m) sets `R_base=Rz(yaw)Ry(pitch)Rx(roll)`, base position in the world, and `params.g=R_baseᵀ params.gWorld`. The dynamics are computed in the base frame; [armS_pressure_geometry.m](armS_pressure_geometry.m) applies the base pose for world-frame animation. Translation of a rigid base changes gravitational potential by a constant, not its gradient. The checked-in runner currently uses three `0.178 m` backbones, `0.013 m` radial offsets, `0.10 kg` per section, a 180° base pitch, `params.gWorld=[0;0;-9.81] m/s²`, and a nonzero abstract generalized input in coordinate 1. Treat those as an **editable example**, not identified hardware values. The separate [simulate_armS_pressure_passive.m](simulate_armS_pressure_passive.m) uses different section lengths, masses, damping, zero input, and a 2 s run.

## 3. Local kinematics already available

For a section at normalized material coordinate `xi∈[0,1]`, [HTM_nume_mex.m](HTM_nume_mex.m) returns `T_i(xi), R_i(xi), P_i(xi)` from `(p_i,xi,L_i_vector,r_i,K_i,A_i)`. [LocalJacob_nume_mex.m](LocalJacob_nume_mex.m) returns `∂P/∂q_i` (3×2), stacked `∂R/∂q_i` (3×6), and position/rotation second derivatives (6×2 and 6×6). These derivatives are **already with respect to the two Pa-valued local coordinates**; applying a length-to-pressure chain rule again would be wrong.

The older helper [pressure2config.m](pressure2config.m) encodes, for `p=[0,p2,p3]`,

```text
alpha = atan2[-sqrt(3)*(p2-p3), 3*(p2+p3)]
beta  = sqrt(3)*A*r*sqrt(p2^2+p2*p3+p3^2)/K
```

Here `alpha` is a bending-plane direction and `beta` behaves as the **section-tip rotation angle** of the exported transform. At `p=[0;2000;1000] Pa`, `Lbody=[0;0.178;0] m`, `r=0.013 m`, `A=pi*(r/2)^2 m²`, and the runner's example `K=0.01443936`, the helper gives `beta=0.547622951`, and `acos((trace(R_tip)-1)/2)=0.547622951` rad. Thus the tested curvature magnitude for a constant-curvature interpretation is `beta/L_i=3.07653343 m⁻¹`. This comparison is at one pose; verify axis/sign conventions and constant strain over the intended range before imposing them in a rod model. `alpha` is undefined at the straight shape; use a two-component bending vector or a strain chart that remains regular there. The helper wraps `beta` to `[-pi,pi]`; do not use that wrapped output as a globally smooth dynamic state.

The older [pressure2length.m](pressure2length.m) expression algebraically simplifies away from its removable zero-pressure `0/0` to

```text
Delta ell_2 = (3*A*r^2/(2*K))*p2,
Delta ell_3 = (3*A*r^2/(2*K))*p3,
Delta ell_1 = -(Delta ell_2+Delta ell_3)   [documented gauge constraint].
```

That helper is **not used** in the three-section core; it is a useful cross-check and possible coordinate bridge. The simplified relation gives zero continuously at `p2=p3=0`, whereas the currently expanded helper should not be called exactly there without fixing its removable singularity. Its physical validity under gravity, external loading, or transient actuation is unverified. Do not identify its `p` with measured chamber pressure during passive motion.

The symbolic sources are in the sibling `Dyn_Pressure/Maple` folder, notably `20261002_3_recursive_T13.mw`, `20261002_3_recursive_extension_T13.mw`, and `20251125_2_pressure_dynamic_model_2DoF_MSF.mw`. The checked-in `*_nume.m` exports are the currently executed local functions. The `_mex` suffix on two `.m` filenames does **not** mean that a compiled three-section MEX exists.

For connected sections, the present backbone model composes

```text
x_i(xi,q) = b_i(q_<i) + R_base,i(q_<i) P_i(xi,q_i),
b_(i+1)  = b_i + R_base,i P_i(1,q_i),
R_base,i+1 = R_base,i R_i(1,q_i).
```

All terms in these physical connections use `Lbody_i`. Full `L(:,i)` is reserved for a future, explicitly defined sensor measurement transform. The present renderer does not establish sensor attachment frames or mounting rotations.

## 4. Existing dynamics and energy

The interpreted core [armS_pressure_core.m](armS_pressure_core.m) returns `M(q)` (6×6), `C(q,qdot)` (6×6), `G(q)` (6×1), and `dM(:,:,h)=∂M/∂q_h` (6×6×6). Section mass is uniform in normalized `xi`: `dm_i=m_i dxi`, equivalently line density `m_i/L_i` kg/m for `s=L_i xi`. For a material point `x_i` and `J_i=∂x_i/∂q`, its modeled energies are

```text
T = 1/2 qdotᵀ M(q) qdot,
M(q) = sum_i m_i ∫_0^1 J_i(xi,q)ᵀ J_i(xi,q) dxi,
V_g(q) = -sum_i m_i ∫_0^1 gᵀ x_i(xi,q) dxi,
G(q) = ∂V_g/∂q.
```

Upstream section rotations contribute to downstream **translational** point velocities. The local-local mass block uses `RᵀR≈I`, so it is not algebraically exact for a truncated polynomial rotation at arbitrary bend. `C` uses the Christoffel convention

```text
C_ab = 1/2 sum_h [M_ab,h + M_ah,b - M_bh,a] qdot_h.
```

Four generated section integrals let the current core avoid repeated point quadrature. For local position `P(xi)` and local pressure derivatives `P_,a`:

| MATLAB export | Moments over `xi∈[0,1]` |
|---|---|
| [integratedPosition_nume.m](integratedPosition_nume.m) | `mu=∫P`, `mu_q`, `mu_qq` |
| [integratedPositionProduct_nume.m](integratedPositionProduct_nume.m) | `S=∫PPᵀ`, `S_q` |
| [integratedPositionDerivativeProduct_nume.m](integratedPositionDerivativeProduct_nume.m) | `F_a=∫P_,a Pᵀ`, `F_q` |
| [integratedJacobianProduct_nume.m](integratedJacobianProduct_nume.m) | `E_ab=∫P_,aᵀP_,b`, `E_q` |

These are integral **moments**, not products of separate integrals: generally `S≠mu muᵀ` and `E≠mu_qᵀmu_q`. The reference recursion is detailed in [recursive_dynamics_pressure_coordinate_handoff.md](recursive_dynamics_pressure_coordinate_handoff.md) and the implemented loops in `armS_pressure_core.m`. A Cosserat implementation need not copy the huge expanded scalar expressions; independent point quadrature is appropriate for the first rod benchmark.

The simulation [armS_pressure_dynamics.m](armS_pressure_dynamics.m) solves

```text
M(q) qdd + C(q,qdot) qdot + G(q) + H q + D qdot = Q.
```

The pressure-coordinate elastic law [pressure_elastic_force.m](pressure_elastic_force.m) is inherited from the existing one-section `G_and_K_pressure_2DoF_MSF` at zero gravity. For section `i`, define `a_i=A_i²r_i²/K_i` and `b_i=9r_i²/(4K_i)`; its 2×2 block is

```text
H_i = a_i [ 3+b_i*k(2,i), 3/2;
            3/2,          3+b_i*k(3,i) ],
U_i = 1/2 q_iᵀ H_i q_i,       F_elastic,i = H_i q_i.
```

`k(1,i)` enters the runner's example choice `K_i=(3/20)k(1,i)L_i r_i²`, not this block directly. `D` is an illustrative 6×6 positive semidefinite pressure-coordinate damping matrix. The energy diagnostic [armS_pressure_energy.m](armS_pressure_energy.m) uses `E=T+sum_i U_i+V_g`. For `Q=0`, a consistent passive model satisfies `dE/dt=-qdotᵀD qdot`; for general `Q`, `dE/dt=qdotᵀQ-qdotᵀD qdot` under the stated model. These are model consistency identities, not experimental validation of `H`, `D`, or `Q`.

## 5. Physical actuation and the critical pressure distinction

The unit of `Q` is `J/Pa=m³` because `δW=Qᵀδq`. It is **not** a pressure in Pa, and the runner's `inputForce` is not a regulator command. A vented arm may change shape while measured chamber gauge pressure remains near zero; a sealed chamber may change pressure as its volume changes. The current `q` traces do not predict either measurement.

For a physically based rod, keep mechanical strain/shape and actual chamber gauge pressure `p_ch` separate. A candidate chamber work model is `δW_p = p_chᵀ δV`, giving `Q_shape=(∂V/∂shape)ᵀp_ch`. If the pressure-equivalent chart is retained, `Q_q=(∂V/∂q)ᵀp_ch` only after defining and validating `V(q)`. The actuator volume/force law, chamber geometry, common-mode pressure, constraints, valve flow, and chamber thermodynamics are **not** in this repository's three-section model. With three physical PMAs but two shape coordinates per section, do not invent a 9×9 invertible mechanical pressure mass matrix; allocate physical chamber pressures to the available shape efforts and enforce feasible pressure limits. Avoid counting the same pneumatic work in both an elastic potential and an applied load.

## 6. What has actually been checked

MATLAB R2025a checks, recorded in [PROJECT_LOG.md](PROJECT_LOG.md), include:

| Check | Evidence and scope |
|---|---|
| Local transform/Jacobian/Hessian and moments | [validate_pressure_local_models.m](validate_pressure_local_models.m): four local poses, centered pressure differences and 24-node quadrature; maximum reported local derivative scaled error `4.529e-6`, quadrature error `2.914e-7`. The exact straight-pressure removable singularity in `HTM_nume_mex.m` was patched by its analytic limit. |
| Three-section M, C, G | [validate_armS_pressure_core.m](validate_armS_pressure_core.m): three poses, 20-node independent global-point quadrature; mass symmetry and Cholesky positive definiteness; `dM` finite differences; `Mdot-2C` skew property; gravity potential-gradient check; sensor-offset invariance. Maximum reported relative `M` quadrature difference `2.298e-8`. |
| Elastic law and interpreted trajectory | [validate_armS_pressure_dynamics.m](validate_armS_pressure_dynamics.m): inherited elastic expression matched at relative error `9.115e-10`; potential-gradient finite difference at `7.337e-8`; 2 s passive trajectory finite and monotonically decreasing in sampled mechanical energy. Its sampled energy-balance residual was `4.758e-5` of the energy drop. |
| World/base geometry | [validate_armS_pressure_geometry.m](validate_armS_pressure_geometry.m): straight and trajectory poses, connected joints, offset invariance, rotation/translation covariance. Maximum reported base-pose covariance error `4.025e-16`. |

These checks cover chosen numerical poses, not a characterized operating envelope. The passive run reached a `13.642 kPa` equivalent coordinate, above the `10 kPa` local validation pose. At the three core poses the largest reported rotation-orthogonality defect was `8.803e-8`; it may grow at larger bends. No experimental shape, force, pressure, or dynamic data have been compared, and no three-section MEX implementation or performance benchmark exists yet.

## 7. Resolve units and constitutive meaning before fitting a Cosserat rod

**Critical dimensional audit:** The coded bend-angle formula `beta=sqrt(3) A r |p|/K` needs `K` in **N·m** if `A` is m², `r` is m, and `p` is Pa. The numerical tip-rotation check above supports interpreting `beta` as an angle, rather than curvature. The comments in local files instead say `N/rad` (dimensionally N if rad is dimensionless). Moreover, the current example `K=(3/20) k L r²` has units **N·m²** if the documented `k` units really are N/m. At least one definition, coefficient, or unit label therefore needs correction or experimental reinterpretation. The current code can be numerically self-consistent while its parameter units are not. **Do not turn this `K` directly into rod bending rigidity `EI` or fit material properties from it until the original derivation and units are settled.**

The same audit applies to `A` (effective fitting area versus chamber area), `k` (PMA axial stiffness versus distributed stiffness), `r` (centerline-to-actuator path), and the inherited `H` and `D`. The existing model assumes a fixed, inextensible section length for dynamics. Measure or identify actual axial extension, shear, torsion, cross-sectional inertia, hybrid rigid features, and any nonuniform line density before expanding the rod's constitutive state. Sensor offset frames and sensor hardware mass must be specified independently of the backbone length.

## 8. Suggested Cosserat formulation and efficiency path

This section is a **proposal**, not existing code. For section `i`, use physical arc length `s∈[0,L_i]`, centerline `x_i(s,t)`, orientation `R_i(s,t)∈SO(3)`, and body strains

```text
v_i = R_iᵀ ∂x_i/∂s,                    [shear/extension]
u_i = vee(R_iᵀ ∂R_i/∂s).             [bend/twist]
```

In spatial force/moment variables `n_i,m_i`, with line density `lambda_i`, distributed force `f_i`, distributed moment `l_i`, and spatial angular momentum per length `h_i`, the continuum balances to implement are `∂_s n_i+f_i=lambda_i ∂_tt x_i` and `∂_s m_i+x_i,s×n_i+l_i=∂_t h_i`. Set `f_i=lambda_i g+f_act/contact` with a consistent gravity frame. Constitutive laws relate the internal force/moment to the strain and strain rate. The clamped base supplies pose/velocity boundary conditions; at each section junction enforce pose, velocity, and wrench transmission; specify free or loaded tip conditions explicitly. For parity with the present reduced model, use `lambda_i=m_i/L_i`, zero added contact/payload, and omit cross-sectional rotary inertia initially. A later physical rod may need nonzero `h_i`, shear, extension, and torsion.

Start by **extracting strains from the existing local geometry** at several poses: differentiate its `x(s)` and `R(s)` with respect to `s` numerically or symbolically, inspect whether `v≈[0;0;1]`, `u_3≈0`, and the two bending components are nearly constant along each section. This tests whether one piecewise-constant-strain (PCS) element per section exactly reproduces the present geometric approximation. Use `beta/L_i` only as an initial bend-magnitude estimate, and determine bend-axis signs from `RᵀR_s`, not from the helper's variable names. Compare reconstructed points and rotations at several `s` values, including straight and asymmetric bends.

For a first **parity** model, retain the same total masses, uniform line density `m_i/L_i`, gravity, two bending degrees per section, omitted rotary inertia, and offset-free physical backbone. Use exact `SO(3)`/`SE(3)` exponentials for rod kinematics and numerical quadrature for kinetic and potential energies. Determine whether the polynomial reference and exact exponential agree within the tested bend range; exact rotations need not match a truncated polynomial at high bend. Include section-to-section pose continuity and a clamped, rotatable base. Exclude sensor offsets from material length and mass, then add sensor measurement poses separately.

After kinematic parity, introduce physically identified rod constitutive laws (for example `n=K_se(v-v0)+D_se vdot`, `m=K_bt(u-u0)+D_bt udot`) and cross-sectional rotary inertia only when corresponding properties are available. A full shear/extension/torsion rod has more modes and cannot be claimed equivalent to the present 6-DOF model by matching six coordinates alone. Add chamber pressure through a derived volume/virtual-work coupling, with pressure dynamics if transient regulator behavior matters. Use the current `H` and `D` only as *numerical comparison baselines* until calibrated.

For efficiency, start with one PCS element per physical section, then refine to 2–4 elements/section and check convergence of tip pose, mass/energy, and trajectories. Assemble sparse/local element contributions and exploit serial recursion. Cache section quadrature nodes and constant material matrices; generate compiled code only after interpreted parity and profiler evidence. Compare the reduced PCS route with a spatial shooting/implicit integration route if distributed shear, twist, or changing curvature are essential. Renda et al.'s [discrete multisection Cosserat model](https://arxiv.org/abs/1702.03660) and Till et al.'s [real-time spatial Cosserat dynamics](https://journals.sagepub.com/doi/10.1177/0278364919842269) are primary methodological starting points, not sources of this arm's parameter values.

## 9. Concrete handoff tasks and acceptance checks

1. **Lock a physical data sheet:** actual `L_i`, mass and distribution, cross-section, actuator paths, chamber volumes, sensor attachment frames, base orientation, pressure range, and whether chambers are vented or sealed. Resolve the `K`/`k` dimensional conflict before material identification.
2. **Build kinematic parity first:** sample `x_i(s),R_i(s)` from the exported model and extract `u_i(s),v_i(s)`; compare an exact PCS reconstruction to the source at straight and asymmetric poses. Report tip and distributed-point errors and rotation orthogonality.
3. **Build interpreted reduced dynamics:** use the same physical mass/gravity assumptions to compare rod and current `M`, gravity, free response, and energy. If the new rod includes extra physics, quantify those contributions separately rather than forcing an exact match.
4. **Derive actuation separately:** use chamber-volume or measured actuator force data to map actual chamber gauge pressure to rod generalized load. Keep valve/pneumatic states distinct from shape coordinates and enforce pressure limits.
5. **Demonstrate efficiency:** record assembly/integration time, memory, element-count convergence, and MATLAB-versus-compiled numerical parity for representative poses and trajectories. The very large existing Maple moment exports previously caused code-generation stalls in the historical model; avoid compiling all expanded expressions as the first step. See [recursive_dynamics_pressure_coordinate_handoff.md](recursive_dynamics_pressure_coordinate_handoff.md).
6. **Validate against measurements:** tip/shape trajectories, chamber pressures, applied inputs, and loads. Internal energy and derivative identities are necessary checks, not evidence that the physical plant is modeled accurately.

**Useful MATLAB entry points:** `validate_pressure_local_interfaces`, `validate_pressure_local_models`, `validate_armS_pressure_core`, `validate_armS_pressure_dynamics`, `validate_armS_pressure_geometry`, and `runDynamicSimulation_pressure_N3`. Use the tests before replacing source functions. Keep [PROJECT_LOG.md](PROJECT_LOG.md) updated with assumptions, parameters, results, and limitations.
