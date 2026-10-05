# Handoff: three-section recursive dynamics and pressure-coordinate extension

Prepared for a separate Codex chat and repository, 5 October 2026. This describes the MATLAB work on `hybrid-soft-arm-dynamics-multisection-pointmass`, then identifies what must be derived anew when the generalized coordinates are pressures, section lengths are a 3×1 vector, and sensors have geometric offsets.

## 1. Start here in the new repository

1. Read the current code and its tests first. The new repository already has its own `*_nume.m` functions for the recursive algorithm, with **uncompacted** symbolic expressions. It also includes **sensor offset extensions**. Treat those local files as authoritative; the original code reviewed below is a historical reference, not a replacement to copy wholesale.
2. Write down the actual configuration variables, their dimension, and the pressure-to-geometry map before changing an equation. “Pressure coordinates” can mean generalized coordinates or actuator inputs; these have different equations of motion.
3. For every section, use its own `L(i)` from a **3×1 length vector**. The old project used one scalar `L` for all sections. Inspect every kinematic, Jacobian, Hessian, integral, gravity, and code-generation call site for this difference. Distinguish physical section lengths from sensor offsets.
4. Preserve the local code’s sensor geometry. Establish whether an offset is a fixed sensor location, a mass-bearing element, or a change to the backbone geometry before including it in kinetic and potential energy.
5. Reproduce the pressure model in interpreted MATLAB and verify analytic derivatives, mass properties, and direct integration before compiling MEX. Avoid the previous long full-build cycle with expanded `S`, `F`, and `E` expressions.

**Historical source repository:** <https://github.com/DulanjanaPerera/hybrid-soft-arm-dynamics-multisection-pointmass.git>. `Ozi` was the point-mass baseline branch; `compare-pointmass-standard` hosted the distributed-mass comparison. The user later built the standard-model MEX through a local Codex session. This handoff does not claim to know that session’s final unshared source diff or performance results.

## 2. Historical model and conventions

- Three serial continuum sections, two independent length-change coordinates per section. The old coordinate vector was

  `q_l = [l_12; l_13; l_22; l_23; l_32; l_33]` (6×1), and `X=[q_l; qdot_l]` (12×1).
- A section’s local kinematic input was `[0, q_l(2i-1), q_l(2i)]`. The first of three actuator-length entries was fixed as a gauge/reference in that formulation.
- `HTM_nume` returned the local homogeneous transform, rotation, and position at material coordinate `xi`, with `xi∈[0,1]`. `LocalJacob_nume` returned position and rotation derivatives through second order. The available symbolic functions used the selected **13th-order Taylor expansion** of the section transform.
- The old scalar `L` was a common section length; the local radial parameter was `r`. The pressure project’s `L=[L_1;L_2;L_3]` changes the interface and possibly the formulas. Pass `L(i)` to section-local functions only after checking the new function contracts.
- The arm is oriented downward in the setup. The stored gravity vector in the reference runs was `[0;0;-9.81]`. The plotted arm uses a `diag([1,-1,-1])` display rotation. Respect the current repo’s base/world transformation and the force-versus-potential sign convention; do not infer a sign from a plot alone.
- Rotational kinetic energy was intentionally omitted in the comparison. **Upstream rotation still contributes to a material point’s translational velocity** and must remain in its Jacobian.
- The historical point-mass model placed mass at a chosen backbone point `cog_xi(i)`, often 0.5 for comparison. That chosen point was not automatically the center of gravity computed by distributed integration. Empirical kinetic-energy beta coefficients were set to one for the validated baseline.
- Sections had fixed total masses `m_i`; distributed integrals used `dm_i=m_i dxi` with the normalized material coordinate. With a nonuniform density, replace this measure consistently in all integrals.

## 3. Original recursive point-mass implementation

Important MATLAB files in the historical repo:

| File | Role |
|---|---|
| `HTM_nume.m`, `HTM_nume_mex.m` | Section transform at `xi` and at the tip |
| `LocalJacob_nume.m`, `LocalJacob_nume_mex.m` | Local first and second transform derivatives |
| `Mi_h.m` | Section mass-matrix derivative with respect to one generalized coordinate |
| `armS_core_N3_mex.m` | Fixed-three-section recursive accumulation of M, C, gravity |
| `christoffelSymbol.m` | Contraction of mass derivatives with generalized velocity |
| `armS_dynamics_N3_entry_mex.m` | ODE RHS, stiffness, damping and generalized input |
| `build_armS_N3_mex.m` | MEX build using a short path outside OneDrive |

For section `i`, the algorithm computed the local transform at `xi=1` and at the chosen mass location `cog_xi(i)`, plus the corresponding local Jacobians and Hessians. It propagated tip position/orientation and their derivatives from preceding sections; it then assembled section mass and gravity contributions. The implementation used 3×3 rotation-derivative blocks, 3×1 position-derivative columns, and section-local derivative arrays sized 2×2×2, 4×4×4, and 6×6×6 for the three sections. It summed their Christoffel contributions into the global 6×6 C matrix.

The important prior corrections were:

1. Bring `HTM_nume` and local Jacobians into alignment with the selected 13th-order expansion and correct derivative indexing.
2. Include each section mass `m_i` in **both** its M contribution and every corresponding derivative slice. The derivative of `m_i M_i` contains `m_i` when section mass is constant.
3. Correct the Christoffel index order for `dM(j,k,h)=∂M(j,k)/∂q_h`:

   `C(j,k) = 0.5 * sum_h [dM(j,k,h)+dM(j,h,k)-dM(k,h,j)] * qdot(h)`.

4. Keep the energy beta coefficients at one for that corrected run. `Mi_h.m` had beta multipliers while the principal M assembly effectively had unit beta; any later nonunit-beta model must be derived from a **single coherent kinetic energy** so M and dM share the factors. The old discussion about missing beta factors is independent of the new pressure-coordinate transformation.

DJ reported a clear simulation improvement after the M correction and excellent behavior after C was fixed. That is a useful consistency signal, not measurement-based proof of the physical dynamics.

## 4. Distributed-mass reference: the analytical integrals

The old project then built a standard distributed-mass model for comparison with the point-mass model, still with translational kinetic energy only. Let `p_i(xi,q_i)` be a point along section `i`. Four generated MATLAB functions supplied the needed section-local moments and derivatives:

| Function | Returned quantities | Definition and dimensions |
|---|---|---|
| `integratedPosition_nume(l,L,r)` | `mu,mu_q,mu_qq` | `mu=∫p dxi` (3×1), `mu_q=∫p_q dxi` (3×2), `mu_qq=∂mu_q/∂q` (3×2×2) |
| `integratedPositionProduct_nume(l,L,r)` | `S,S_q` | `S=∫ppᵀ dxi` (3×3), `S_q=∂S/∂q` (3×3×2) |
| `integratedPositionDerivativeProduct_nume(l,L,r)` | `F,F_q` | `F(:,:,a)=∫p_{,a}pᵀ dxi` (3×3×2), `F_q(:,:,a,b)=∂F(:,:,a)/∂q_b` (3×3×2×2) |
| `integratedJacobianProduct_nume(l,L,r)` | `E,E_q` | `E=∫p_qᵀp_q dxi` (2×2), `E_q=∂E/∂q` (2×2×2) |

Every integral above is over `xi=0…1`. `a=1` differentiates with respect to `l(2)`; `a=2` differentiates with respect to `l(3)`. Index `b` has the same order. Crucially, `F_a` is `∫p_{,a}pᵀ`, **not** `∫pp_{,a}ᵀ`. The identity `S_{,a}=F_a+F_aᵀ` is a useful implementation check.

The CoG paper’s partial-trace identities may be coded as ordinary MATLAB block products and traces. The identities rearrange matrix products; they do not allow `∫(abᵀ)` to become `(∫a)(∫bᵀ)`. In particular `S≠mu muᵀ` and `E≠mu_qᵀ mu_q` in general. One may differentiate an analytical integral after integration when the integration bounds and mass density are parameter independent, or integrate the full differentiated expression; product rules remain essential.

### Why these quantities appear

Let `P` and `R` be the accumulated position and orientation of the base of section `i` in a common frame. A material point has position

`x_i(q,xi)=P(q_up)+R(q_up) p_i(q_i,xi)`.

For upstream coordinate `j`, define `A_j=P_,j`, `B_j=R_,j`; for local coordinate `a`:

`J_j=A_j+B_j p_i`, `J_a=R p_{i,a}`.

The translational kinetic energy for section `i` is

`T_i = (m_i/2) ∫_0^1 ||J_i(q,xi) qdot||² dxi`.

The historical standard core assembled the following section mass blocks (the section contributes to all preceding coordinates):

```text
M_jk / m_i = A_jᵀ A_k + A_jᵀ B_k mu + A_kᵀ B_j mu
              + tr(B_jᵀ B_k S)
M_ja / m_i = A_jᵀ R mu_q(:,a) + tr(B_jᵀ R F(:,:,a))
M_aj       = M_ja
M_ab / m_i = E(a,b)
```

The local-local block uses the ideal rotation identity `RᵀR=I`. With the truncated Taylor transform, `R` is only approximately orthogonal; the assembled M and directly integrated global point Jacobians agreed very closely at tested configurations, but do not assert exact identity at arbitrarily large bending. If the new pressure implementation has substantial orientation truncation error, derive an exact local-local block including `RᵀR`, use a sufficiently accurate transform, or quantify the approximation.

`armS_standard_core.m` propagated global `P_,j`, `R_,j`, and second derivatives through

`P_next=P+R p_tip`, `R_next=R R_tip`,

using the product rule at each section. It analytically differentiated the block formulas (including `mu_q`, `mu_qq`, `S_q`, `F_q`, `E_q`) to obtain `dM(:,:,h)`. Then it used `christoffelSymbol` for C. The gravity term used the same sign/convention as the point-mass model:

`G_j += m_i (A_j+B_j mu)ᵀ g`, `G_a += m_i (R mu_q(:,a))ᵀ g`.

Do not copy this gravity sign mechanically to the new repo: first confirm its frame and whether its `G` denotes physical potential gradient or the signed force convention used by its ODE RHS.

The old RHS kept the same stiffness, damping, and input law as the corrected point-mass baseline:

`qddot_l = M_l \ [tau_l - (C_l+D_l) qdot_l - G_l - K_l(q_l) q_l]`.

For the new pressure-coordinate model, these length-coordinate stiffness and damping terms require transformation or rederivation from an energy/dissipation law. A pressure command is not automatically a generalized pressure-coordinate force.

## 5. Tests and observed results

For the standard model, `validate_armS_standard.m` checked M symmetry, Cholesky positive definiteness, centered finite differences of each `dM(:,:,h)`, skew symmetry of `Mdot-2C`, and agreement of M and gravity with independent 16-node Gauss–Legendre integration of global material-point Jacobians. DJ provided:

```text
Pose 1: symmetry 0.000e+00, dM 0.000e+00, skew 0.000e+00,
        integral M 1.817e-15, G 0.000e+00
Pose 2: symmetry 1.893e-16, dM 2.151e-10, skew 1.672e-16,
        integral M 1.755e-15, G 1.574e-15
Pose 3: symmetry 8.891e-17, dM 5.462e-11, skew 8.213e-17,
        integral M 1.112e-12, G 1.550e-15
All standard-model checks passed at the three test poses.
```

The interpreted standard simulation finished with finite 301×12 states over five seconds and decaying oscillations. Against the corrected point-mass run using the same saved parameters, initial state, and time grid, the per-coordinate length RMSEs were approximately 1.108, 0.764, and 0.403 mm for sections 1–3; the maximum absolute coordinate difference was approximately 2.520 mm. These measure model disagreement, not error against experimental truth. The observed interpreted integration time was 15.786 seconds; point-mass timing involved compiled code and is not a fair like-for-like benchmark.

## 6. The MATLAB Coder bottleneck and what resolved it

The expanded Maple exports for `S`, `F`, and `E` and their derivatives were large. MATLAB Coder spent many minutes after `### Compiling function(s) ...` without producing source files. MATLAB sometimes did not respond to Ctrl+C and had to be ended from Task Manager. An earlier first build also failed because loops written as `for i=old` had a variable-sized index vector; explicit ranges `for i=1:2*(n-1)` worked for code generation.

To isolate the bottleneck, source-only tests used `cfg=coder.config('mex')`, `cfg.GenCodeOnly=true`, report generation off, `cfg.InlineBetweenUserFunctions='Never'`, and verbose `codegen('-v',...)` with a short `-d` path under `C:\MATLAB_build`. The test of the original `mu` generated source in about 35.09 s. The original `F` and `S` stalled sufficiently long that they were interrupted.

The remedy for S/F/E was to parse their original scalar expressions, use exact rational constants during algebraic rewriting, and eliminate common subexpressions, producing separate functions:

```matlab
integratedPositionProduct_compact(l,L,r)
integratedPositionDerivativeProduct_compact(l,L,r)
integratedJacobianProduct_compact(l,L,r)
```

They represent the **same analytical integral expressions and derivative tensors**, not numerical integration or a different Taylor expansion. Keep the original `_nume.m` functions as references. The user's MATLAB check reported compact `F` and `F_q` matching originals over 35 poses with maximum relative errors `1.857e-15` and `2.430e-15`. Independent assistant-side checks of compact S/S_q and E/E_q covered 70 geometry/pose combinations with relative errors below about `5e-15`. Source-generation times reported by DJ:

| Source-only function | Time |
|---|---:|
| Original `integratedPosition_nume` | 35.091 s |
| Compact F/F_q | 26.038 s |
| Compact S/S_q | 18.877 s |
| Compact E/E_q | 4.329 s |

Source-only generation is not a full C++ compile or a numerical MEX parity test. The user subsequently said local Codex built the standard model successfully; inspect that repo’s actual build files and logs for details.

**New pressure repository implication:** its already-created, uncompacted `*_nume.m` files may again stall MATLAB Coder. Their expressions may differ from the old scalar-length ones because pressure variables, per-section `L_i`, and sensor offsets change geometry and derivatives. **Do not copy the old compact S/F/E files as if they are equivalent to the new exports.** Validate each new original function against direct point sampling/quadrature and derivative identities, then compact those exact new scalar expressions, test compact-versus-original at representative physical configurations, and source-generate one function at a time before an integrated build. Log the exact source version used to build each MEX.

## 7. Pressure coordinates: mathematical decision before coding

Let the new generalized pressure coordinates be `pi` (use another name in MATLAB to avoid shadowing the mathematical constant `pi`), and let the old geometric coordinates be `ell`. First identify the new project’s actual mapping

`ell = f(p)` with `p∈R^m`, `ell∈R^6`.

Define `J_f = ∂ell/∂p` (6×m) and `H_f(:,:,a) = ∂J_f/∂p_a`. If the new `_nume` kinematics are **already differentiated directly with respect to p**, use those derivatives and do not apply this transform twice. If they still use length coordinates internally, use the chain rule:

```text
ell_dot  = J_f p_dot
ell_ddot = J_f p_ddot + Jdot_f p_dot
Jdot_f   = sum_a H_f(:,:,a) p_dot(a)
```

For a material point, `J_x,p = J_x,ell J_f`; its second derivatives need both geometric Hessians and mapping Hessians. In components,

`∂²x/∂p_a∂p_b = Σ_jk x_,ell_j ell_k J_f(j,a) J_f(k,b) + Σ_j x_,ell_j ∂²ell_j/∂p_a∂p_b`.

For **a holonomic coordinate reparameterization with an invertible, square mapping** (here m=6 and nonsingular J_f), transform kinetic and generalized forces consistently:

```text
M_p     = J_fᵀ M_ell J_f
G_p     = J_fᵀ G_ell                       [with the same sign convention]
Q_p     = J_fᵀ Q_ell                       [work conjugacy]
h_p     = J_fᵀ [ M_ell (Jdot_f p_dot)
               + C_ell(ell,ell_dot) ell_dot ]
```

The velocity-quadratic vector `h_p` can be represented as `C_p(p,p_dot)p_dot`; the matrix `C_p` is not unique, but it must reproduce that vector and satisfy the standard mass-derivative consistency checks if formed by Christoffel symbols. Elastic energy, damping, and actuator work must also be transformed from their defining energies/work laws. For example, for a length-coordinate viscous generalized force `D_ell ell_dot`, the transformed force is `J_fᵀ D_ell J_f p_dot`; for potential `V(ell(p))`, its gradient is `J_fᵀ ∇_ell V`.

**Dimensionality warning:** if the new system treats all three physical PMA pressures in each of three sections as independent generalized coordinates, then m=9 while the old backbone model has only six geometric coordinates. The induced `M_p=J_fᵀ M_ell J_f` has rank at most six and cannot be inverted as a 9×9 mass matrix. The three redundant directions require constraints, a reduced pressure chart, or additional physical states and energies (for example pneumatic dynamics and fluid/compliance energy). A 3×1 vector `L` denotes section lengths and does **not** establish whether the pressure-coordinate dimension is 6 or 9. Inspect the actual new repository for that dimension. Never hide this structural singularity with arbitrary diagonal regularization.

If pressures are **inputs** to a length-coordinate mechanical model rather than coordinates, keep the mechanical coordinates as lengths/shape variables and derive the actuator generalized forces `Q_l(p,ell)` plus pressure dynamics if needed. There is no automatic `M_p` for a pressure command vector. Clarify whether the model’s pressure signals are commanded, regulated measured pressures, chamber states, or a reduced set of shape pressures.

## 8. Vector section lengths and sensor offsets

### `L` as a 3×1 parameter

In section loop `i=1:3`, local transforms and integrals should generally receive the section’s own scalar `L(i)`, if their signature means physical length. Examples to audit:

```matlab
HTM_nume(localConfiguration, xi, L(i), r)
LocalJacob_nume(localConfiguration, xi, L(i), r)
integratedPosition*_nume(localConfiguration, L(i), r)
```

This is an **interface illustration**, not a claim that the new repo’s signatures match the old ones. If a new function intentionally takes the full vector to describe cross-section geometry, follow its documented contract. Update code generation argument examples to use a fixed-size 3×1 `L` input, and test nonidentical values, such as `[0.278;0.260;0.290]`, so accidental use of `L(1)` for all sections is visible. If `L` is constant in time, its derivatives with respect to pressure are zero; if physical extension makes `L=L(p,t)`, include its pressure/time derivatives and the appropriate moving-mass or variable-length modeling assumptions.

### Sensor geometry

A sensor offset needs an explicit frame and attachment point. For a sensor with a fixed offset `d_s` in the section frame at material point `xi_s`, a typical point model is

`x_s = P_base + R_base [ p_i(xi_s,p) + R_i(xi_s,p) d_s ]`.

Check whether the local sensor frame includes an additional fixed mounting rotation; if so apply it in the correct order. For fixed local offset `d_s`, the local first derivative includes

`(p_i)_,a + (R_i)_,a d_s`,

and the second derivative includes `(p_i)_,ab + (R_i)_,ab d_s` (plus other terms when the sensor attachment or offset varies with p). Compose these with upstream transform derivatives. Sensor outputs may also be full poses; differentiate the orientation output consistently when needed for measurement Jacobians.

The sensor’s geometric offset does **not** automatically change the material backbone integration or the section’s physical length. If sensor hardware has non-negligible mass, add its translational kinetic/potential energy explicitly and decide whether its rotational kinetic energy is included in the new model. This is separate from the earlier decision to omit rotation energy for the arm. If the offset describes an actual extended arm segment, document that mechanical extension and its mass distribution, then rederive its integrals and section-tip transform. Never fold an arbitrary sensing offset into `L(i)` silently.

A sensor output used only for measurement should not modify M, C, or G. A physical tip payload, extra spacer, or off-axis mass does modify these terms and requires explicit mass and inertia assumptions.

## 9. Validation matrix for the new project

Use independent numerical point sampling of the new pressure-coordinate geometry as the reference, not merely code-to-code comparisons of two functions with a shared error.

1. **Interface/units:** Verify pressure units (Pa versus bar), `L` order and dimensions, sensor-offset frames and units, coordinate/state ordering, total section masses and density measure.
2. **Local functions:** Finite-difference or complex-step check transforms, local position/rotation Jacobians and Hessians with respect to actual pressure coordinates. Validate `mu`, `S`, `F`, `E` by quadrature of the material-point model and check `S_,a=F_a+F_aᵀ`. For any pressure-dependent density or moving integration bound, include derivatives of those terms.
3. **Recursive composition:** Compare each section-tip and selected material/sensor point global position and Jacobian with finite differences of a direct product of transforms. Verify second derivatives, including sensor offsets and pressure-to-length Hessians where relevant.
4. **Dynamics:** Compare recursively assembled M and G with direct quadrature of global point Jacobians; verify M symmetry, finite values, and positive definiteness **only if the chosen coordinates are independent and the modeled kinetic energy can support it**. For redundant pressure coordinates, expect semidefiniteness and address the model dimension explicitly.
5. **Velocity terms:** Compare `dM(:,:,h)` with finite differences and `C p_dot` with a direct Lagrangian/Christoffel construction. Check `Mdot-2C` skew symmetry for a consistent Christoffel convention. Keep gravitational, elastic, damping and input signs documented separately.
6. **Geometry cases:** Include straight and asymmetric bending, unequal `L(i)`, distinct pressures by section, zero and nonzero sensor offsets, and a configuration near the intended operating limit. Check that zero offset reduces to the expected baseline measurement geometry.
7. **Code generation:** Run source-only codegen for each new large analytic integral, compact *that exact function* if required, compare originals to compact versions, then compile a single fixed-interface RHS. Keep physical pressures, section lengths, masses, gravity and sensor offsets as runtime inputs when tuning is intended. Run compiled-versus-MATLAB RHS parity at nontrivial values and trajectory parity before benchmarking.
8. **Simulation comparison:** Use identical physical inputs, initial mechanical state, gravity, mass, stiffness, damping, and solver settings for any length- versus pressure-coordinate comparison. If state dimensions differ, compare physical tip/sensor poses and energy rather than coordinate vectors directly.

## 10. Suggested instruction to the new Codex chat

> Read this handoff and inspect the pressure-dynamics repository, including all existing `*_nume.m` exports and sensor-offset code. First identify the exact generalized coordinate dimension, pressure-to-geometry map, meaning of the 3×1 `L`, sensor-offset frames, mass assumptions, and what physical quantity is used as input. Show the resulting equations and the mapping from each existing MATLAB file to each needed term. Then implement the recursive translational dynamics and its validation in small reviewable steps. Preserve original symbolic exports, compact only the new expressions proven equivalent to those exports, and isolate code-generation bottlenecks before a full MEX build. Record any assumptions that the repository cannot establish.

## 11. Source and confidence boundary

The historical conclusions above come from the MATLAB files inspected during this conversation, the attached CoG and spatial-dynamics references used in the modeling discussion, DJ’s MATLAB validation output, and the saved comparison data. **The new pressure-coordinate repository was not supplied or inspected here.** Its generalized coordinates, exact actuation mapping, sensor offsets and mass distribution must be established from that repository before implementing the formulas. Statements in sections 7–9 about the new model are conditional derivations and implementation guidance, not claims that a working pressure-coordinate recursive model already exists.
