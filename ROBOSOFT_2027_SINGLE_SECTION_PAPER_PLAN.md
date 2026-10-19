# RoboSoft 2027: single-section dynamic-model paper plan

Status: working plan, October 19, 2026. This plan concerns **one hybrid continuum-arm section actuated by three PMAs**, with two bending degrees of freedom. No experimental results are claimed here.

## Submission constraints

- [RoboSoft 2027 contributed-paper deadline](https://robosoft27.com/call-for-papers/): October 31, 2026, 11:59 PM Pacific Time. Internal target: a complete, reviewed manuscript by **Sunday, October 25**; October 26–30 are reserved for corrections and submission checks.
- Six US-letter IEEE-format pages **including references**; up to two paid extra pages are allowed. Draft for six pages. The call encourages a supplementary video (MP4, at most 3 minutes and 20 MB).
- Double-anonymous review: anonymize manuscript, figures, videos, file metadata, and repository links. Follow the call's AI-content disclosure policy in the submitted version as applicable, while preserving review anonymity.

## Research question and claim

**Question:** Can a compact, physically structured two-DOF dynamic model predict the measured tip motion of this three-PMA module from commanded regulator pressures, and improve tracking over simpler command-to-shape control?

**Candidate contribution, conditional on results:** An identified **command-to-motion** dynamic model built on the existing single-section equations, with held-out prediction and a controlled tracking experiment on the physical module. Report the tested command/motion envelope and **effective** fitted stiffness/damping. Chamber pressure is not recorded, so this is not a validated pressure-to-mechanical-effort model and the fitted parameters should not be labeled independently measured physical material constants. Do not claim the first recursive model, the first pressure-coordinate model, or a complete regulator/valve-flow model.

**Available recordings (researcher, October 19):** three **commanded pressure** setpoints, tip `X,Y,Z`, and measured `theta,phi`. There is no chamber-pressure measurement. Confirm synchronized timestamps and sample rates before identification. Tip XYZ and `theta,phi` describe the same two-DOF configuration and should not be treated as five independent mechanical degrees of freedom; use angle trajectories as primary fit outputs and XYZ as an independent geometry/prediction check when those measurements are independently obtained.

**Logging handoff:** The researcher is setting up MATLAB/Simulink logging before sharing a trial. Use [ROBOSOFT_2027_DATA_LOG_SPEC.md](ROBOSOFT_2027_DATA_LOG_SPEC.md) for command/pose columns, units, timestamps, metadata, commissioning steps, and the [save_robosoft_trial.m](save_robosoft_trial.m) table export helper. First obtain one short commissioning trial and inspect it before running chirps.

**Naming contract:** `p_cmd ∈ R^3` is the commanded pressure setpoint for the three pressure regulators; `p_ch ∈ R^3` is actual, **unmeasured** chamber pressure; `y=[theta,phi,X,Y,Z]` is the recorded tip pose; `q_eq ∈ R^2` is the existing pressure-equivalent *shape coordinate*; `Q ∈ R^2` is generalized effort. The existing `arm_dynamics_pressure_2DoF_MSF.m` uses `tau=Q`, not `p_cmd` or `p_ch`. The zero-sum relation used in the shape chart is not a physical constraint on three nonnegative PMA commands.

**Model bridge to establish:**

```text
recorded p_cmd  -->  unknown regulator + PMA response  -->  effective Q
                                                     -->  M(q_eq) q_eq_dd + C q_eq_d
                                                          + G + F_el + D q_eq_d = Q
                                                     -->  predicted theta, phi, XYZ
```

Fit a parsimonious **effective command-to-effort** map, or a command-to-shape model if the physical bridge cannot be justified. A single empirical lag may be included only if it improves held-out prediction; without `p_ch`, its time constant cannot be assigned uniquely to the regulator rather than the PMA or mechanical dissipation. Avoid fitting an unconstrained scale into both the actuation gain and stiffness without an identifiability check. Forward-simulate `q_eq`, then map it to predicted `theta,phi,XYZ`; do not convert measured angles to `q_eq` once using a fixed `K` and then independently optimize `K`. The current definition of `K` and its units must be reconciled before reporting it as physical bending rigidity; otherwise report identified **effective** parameters and their coordinate convention.

## Six-page manuscript layout

| Space | Section | Required content and evidence |
| --- | --- | --- |
| Page 1 | Title, abstract, introduction | Robot photo/diagram; problem; precise contribution; one-sentence numerical result placeholders to fill only after experiments; related-work distinction. |
| Page 2 | Module and model | Three-PMA hardware and two bending coordinates; recorded commands and tip pose; equation of motion; `p_cmd`, unmeasured `p_ch`, `q_eq`, `Q` distinction; assumptions and operating range. |
| Page 3 | Command coupling and identification | Effective command-input map, excitation design, parameter vector, bounded least-squares objective, training/validation split, identifiability statement. |
| Page 4 | Dynamic-model validation | Held-out measured-vs-predicted `theta,phi` and XYZ; static-map and simple dynamic baselines; commanded-pressure traces; error table and repeat variability. |
| Page 5 | Control demonstration | Model-based command feedforward plus measured-pose feedback, command allocation and saturation handling; compare with direct static command mapping and tuned feedback-only control on the same unseen references. |
| Page 6 | Discussion, limitations, conclusion, references | State where the model succeeds/fails, computational cost, uncertainty from unmeasured chamber pressure and hysteresis; concise conclusion and references. Reserve enough space for references from the first draft. |

**Provisional title:** “Identification and Experimental Validation of Command-to-Motion Dynamics for a Hybrid Soft Continuum Module.” Revise after results; include “control” in the title only if the hardware control comparison is convincing.

**Figure/table budget:** (1) module and available sensing schematic; (2) model/input notation and command-to-angle response; (3) held-out two-axis trajectory overlay, XYZ error, and three commanded-pressure traces; (4) control tracking overlay; (5) compact parameter/error/timing table. Combine panels to stay within six pages. Supplementary video should show real hardware and synchronized overlays if available; an animation alone is not experimental validation.

## Minimum executable experiment set

1. **Instrumentation check:** Synchronize timestamps for the three commanded pressures, `theta,phi`, and tip XYZ. Record sampling rates, angular conventions, pose calibration, command limits, and a zero-command baseline. Save raw data and metadata before processing. Document that chamber pressure is unavailable.
2. **System response check:** Use modest command steps and slow sweeps to determine safe amplitude, visible lag, settling time, and inflation/venting asymmetry from **command-to-pose** measurements. Choose chirp frequencies within the observable pose bandwidth; do not label the inferred lag as a measured regulator time constant.
3. **Identification:** Acquire slow/static holds at several feasible command combinations. Excite two independent differential-command directions around a nonnegative bias; test common command as a separate condition because it may alter length or apparent stiffness. Fit a compact shared parameter set with `lsqnonlin` or equivalent bounded optimization on forward-simulated angle trajectories. Use the recorded command as input; compute predicted XYZ from the model as an additional validation output. Hold measured module mass and geometry fixed where credible; estimate only parameters the data can support. Repeat at least one trial to assess variability. A zero-command release/free-decay trial can help constrain effective mechanical damping, but without chamber-pressure sensing it does not prove that the PMAs are vented immediately.
4. **Held-out validation:** Freeze parameters; predict new waveform(s), amplitude(s), and initial pose(s). Include at least one combined-axis trial and a return/venting phase. Report RMSE or normalized RMSE in `theta,phi` and tip XYZ, phase/timing error, and visible failure cases. Compare with a quasi-static **command-to-shape** map and a low-order dynamic baseline using the same training data.
5. **Control:** If the identified model passes validation, run matched reference trajectories with (a) direct static command mapping, (b) feedback-only pose control, and (c) model-based command feedforward plus the same feedback. Use measured pose feedback and respect command bounds. Report tracking error, settling/phase, peak commanded pressure, saturation, and repeated-trial variation. Do not claim physical pressure tracking or exact mechanical feedback linearization without chamber-pressure and actuation-map evidence.

**Data split:** Name training trials before fitting; keep validation and control reference trials untouched until parameter freeze. Keep raw data, processing scripts, fit settings, and plots reproducible. A simulation video from `runDynamicSimulation_pressure_2DoF_MSF.m` is illustrative only.

## Schedule to the Sunday manuscript freeze

| Date (Chicago time) | Deliverable and decision |
| --- | --- |
| Mon Oct 19 | Freeze paper question, coordinate/input notation, available experimental channels, six-page outline, and figure list. Confirm command and pose can be timestamped together. |
| Tue Oct 20 | Make command/pose logging reliable using the data-log specification; inspect one commissioning trial, pose/clock calibration, and basic command-to-pose response. Draft introduction, hardware, and model sections with placeholder-free methods. |
| Wed Oct 21 | Complete identification trials in two bending directions plus static holds; archive raw data. **Gate 1:** if synchronized command and pose data are not usable, narrow the claim immediately. |
| Thu Oct 22 | Fit parameters, check residuals and parameter correlation, freeze training data and parameter set. Draft identification section. |
| Fri Oct 23 | Run held-out dynamic validation and baseline comparisons; generate main validation figure/table. **Gate 2:** if prediction does not beat simple baselines, revise the model or paper claim before control work. |
| Sat Oct 24 | Run matched control comparisons if the model is ready; otherwise devote the day to additional validation and a modeling-only paper. Produce control figure/video if supported. |
| Sun Oct 25 | Complete six-page manuscript, references, abstract with measured numbers, internal technical review, anonymization pass, and PDF/figure checks. Freeze a reviewable version. |
| Oct 26–30 | Address review findings, rerun only experiments needed to resolve concrete claims, check PaperPlaza compliance, upload early. Official deadline: Oct 31, 11:59 PM Pacific. |

## Responsibilities and current unknowns

- **Researcher/hardware team:** Physical excitation and safety limits, tip-pose calibration, synchronized command/pose recordings, controller deployment, and experimental sign-off.
- **Codex/code and writing support:** Audit the single-section equations and units; implement reproducible import, identification, simulation, comparison, and plotting scripts once the data format is known; draft and revise the paper from verified results; maintain `PROJECT_LOG.md` with decisions and measured evidence.
- **Confirmed:** Only commanded pressure, tip XYZ, and `theta,phi` can be recorded; chamber-pressure traces are unavailable. Acquisition is in MATLAB or Simulink. **Still to confirm:** Logging rate/timestamp synchronization, command units and limits, pose measurement method, experiment time, and who can deploy the physical controller. These determine whether the control result is achievable by October 24.

## Paper integrity gates

- Report estimated **effective** stiffness/damping and command-to-motion behavior. Without actual pressure or force sensing, avoid claiming that the fit separately identifies intrinsic mechanical damping, material stiffness, or regulator dynamics. Do not recycle simulation parameters as identified measurements.
- State the pressure-coordinate chart and its dependence on `K`; derive outputs from independent measured shape. Resolve the current `K` dimensional inconsistency before claiming a physical modulus or rigidity.
- State limitations from hysteresis, common pressure, extension, and unmodeled regulator dynamics. Keep all numerical claims traceable to frozen scripts and raw trials.
- If the weekend evidence supports only identification and prediction, submit a modeling/validation paper only if that contribution is sufficiently distinct from prior work; otherwise use the February 2027 non-proceedings extended abstract or a later venue.
