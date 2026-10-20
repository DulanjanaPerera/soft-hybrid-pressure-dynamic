# Important decisions

## 2026-10-05 — Coordinate and mass contracts
Use two pressure-equivalent shape coordinates per section and pressure-local derivatives directly, rather than importing length-coordinate derivatives without transformation. The local leading zero is a reference slot. Only L(2,i) carries backbone mass and contributes physical inter-section length; sensor offsets are massless. Distributed mass uses dm_i=m_i dxi. Rotary inertia is excluded. These choices determine assembly and geometry; future hardware mass requires explicit additions.

## 2026-10-05 — Mechanical effort and chamber pressure are separate
Passive vented motion can change shape while gauge pressure remains near zero. Therefore q is a shape chart, Q is J/Pa, and chamber pressure is a separate physical quantity. A physical pressure-to-effort relation requires derivation/calibration; passing Q directly to a regulator is invalid. Direct inverse pressure mapping remains an experimental control baseline.

## 2026-10-05 — Base/world convention
Use Rbase=Rz(yaw)Ry(pitch)Rx(roll) for rendering and g=Rbase' gWorld for dynamics. Translate the base only in world geometry/potential reference. Gravity G in the pressure core is a potential gradient.

## 2026-10-07 — Compact MEX path
Expanded generated functions stalled Coder analysis. Use compact pressure-local Taylor derivatives and 16-node quadrature for a separate compiled recursive assembly, retaining symbolic MATLAB as reference. Implication: keep both assemblies synchronized and recheck parity over any enlarged range.

## 2026-10-14 — Pressure regulators define the hardware boundary
The researcher confirmed 0–10 V commands drive pressure regulators. Begin with observed command response; explicit raw flow-valve modeling is optional if evidence requires it. Do not assume a first-order lag or independently measured regulator dynamics.

## 2026-10-19 — Immediate paper scope changed to one module
Supersedes the earlier three-section experimental/control paper recommendation. Study one three-PMA module, identify dynamics with informative excitation, validate held-out prediction, and demonstrate control conditionally. Rationale: narrow the experimental contribution to an executable scope. Three-section numerical work remains historical implementation evidence.

## 2026-10-19 — Available channels narrow the claim
Supersedes earlier plans requiring measured chamber pressures. Available recordings are three commanded pressure setpoints, tip XYZ, and theta,phi; no actual chamber-pressure measurements. Fit a command-to-motion model with effective parameters. Regulator/PMA/mechanical lag cannot be separated uniquely from these outputs alone. Fit predicted angles directly; if K varies, recompute its coordinate chart within fitting or use independent mechanical coordinates. Separate trials for training and held-out validation.

## 2026-10-19 — Manuscript delivery
The historical log records that the sibling pressure-control-paper repository uses main for Overleaf, with main.tex as the manuscript and IEEEtran.cls retained. Future paper work follows that destination; no paper changes were made here.

## 2026-10-20 — Documentation contract
Adopt the seven-file project_docs structure requested in GPT_prompt_Log.txt. PROJECT_STATE.md is the primary current handoff; PROJECT_LOG.md is append-only history; specialized files hold decisions, useful results, active issues, architecture, and reusable notes. Keep original handoff reports and experiment specifications as linked detailed sources.
