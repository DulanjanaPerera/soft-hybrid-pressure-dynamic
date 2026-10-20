# Technical notes

Reusable equations and interpretation are collected here. Detailed source derivations and historical parameter tables remain in the [Cosserat handoff](../cosserat_rod_pressure_model_handoff.md) and dated [history](PROJECT_LOG.md).

## Shape inverse and smooth references
In the legacy pressure2config convention, theta is plane direction and phi is bend magnitude:
theta=atan2(-sqrt(3)(q2-q3),3(q2+q3)).
Let s=K/(A r), u=phi cos(theta), v=phi sin(theta). Then
q2=s(u/3-v/sqrt(3)), q3=s(u/3+v/sqrt(3)).
Use 0<=phi<pi for the supported inverse branch. Bend components remain defined at straight shape, where plane angle is undefined. config2pressure_reference_N3 analytically differentiates smooth references; its output is a shape reference, not physical pressure allocation.

## Length chart and pneumatic virtual work
For one section, c=3Ar²/(2K):
ell=c[-q2-q3;q2;q3], J_ell,q=c[-1,-1;1,0;0,1].
The straight singularity in the legacy helper is removable algebraically. This is a coordinate relation, not a causal chamber-pressure law.

If chamber volume V(q) is known and deltaW=p_ch' deltaV applies with consistent gauge pressures and PMA signs, Q=(dV/dq)' p_ch. Under the additional unverified constant-volume-area assumption deltaV=A_V deltaell:
Q=A_V c [p_ch,2-p_ch,1;p_ch,3-p_ch,1].
Common pressure is null only under this simplified assumption; actual PMA extension and stiffness can depend on it. Avoid double-counting pressure work in elastic and actuation terms.

## Elastic and energy conventions
For section i, a=A_i² r_i²/K_i, b=9r_i²/(4K_i):
H_i=a*[3+b*k(2,i),3/2;3/2,3+b*k(3,i)],
F_el,i=H_i q_i, U_i=0.5 q_i' H_i q_i.
This is inherited from the single-section expression, not a calibrated pneumatic constitutive law. With G=grad(V), total E=0.5 dq' M dq+U+V has passive dissipation rate -dq' D dq under the ideal modeled assumptions.

## Ideal computed-effort controller
Let e=q_d-q, zdot=e, h=C dq+G+F_el+D dq.
Q*=M(q)(qdd_d+K_D edot+K_P e+K_I z)+h.
Under exact model, state availability, nonsingular M, and instant feasible application of Q*, the error satisfies edd+K_D edot+K_P e+K_I z=0. For independent diagonal axes, positive gains and K_D*K_P>K_I componentwise give the cubic Routh-Hurwitz condition. Arbitrary positive matrix gains do not suffice. This is not a validated regulator-level controller.

## Identifiability
Illustration: I theta_dd+d theta_d+k theta=a DeltaP. Static slope identifies a/k. Known I and informative dynamics can in principle separate k,d,a; freeing I permits a shared scale ambiguity. In the current setup p_ch is not recorded, so use p_cmd with an explicitly effective input model. A fitted lag cannot uniquely identify regulator lag. Since the q chart contains K, fitting K after fixing a shape-to-q conversion is inconsistent; recompute the chart or fit measured mechanical coordinates.

## Units and frames
The angle formula sqrt(3) A r sqrt(q2²+q2 q3+q3²)/K requires K in N·m for q in Pa. The inherited K=(3/20)k L r² has N·m² when k is N/m. This unresolved discrepancy is retained rather than silently repaired.
For base Euler angles, Rbase=Rz(yaw)Ry(pitch)Rx(roll); g_base=Rbase' g_world. Sensor offsets are excluded from material propagation/mass. Historical length-model gravity sign differs from the pressure potential-gradient convention.

## Numerical implementation
Compact MEX position uses the exported transform's polynomial through power 12; rotation uses a 30-term exponential series. Sixteen Gauss-Legendre nodes integrate position-product polynomials of degree at most 26 exactly in exact arithmetic. Floating-point evaluation and separately generated symbolic integrals still differ, especially at larger bends. At local [0,10000,0] Pa, one symbolic first-moment component differed from transform quadrature by about 7.4e-8 m. Selected-pose parity is not an operating-envelope error bound.

## Measurement and experiment interpretation
Fit forward-simulated angle trajectories, preserving wrap conventions and straight-pose ambiguity. Use XYZ as an additional geometry check; outputs from the same NDI pipeline are correlated. Keep raw streams and device timestamps, document alignment/filtering, and split by whole trial. No measured chamber pressures means mechanical and regulator/PMA effects cannot be isolated uniquely. The detailed [logging specification](../ROBOSOFT_2027_DATA_LOG_SPEC.md) remains the acquisition schema authority.
