# RoboSoft 2027: single-section dynamic-model paper plan

Status: working plan, October 19, 2026. This plan concerns **one hybrid continuum-arm section actuated by three PMAs**, with two bending degrees of freedom. No experimental results are claimed here.

## Submission constraints

- [RoboSoft 2027 contributed-paper deadline](https://robosoft27.com/call-for-papers/): October 31, 2026, 11:59 PM Pacific Time. Internal target: a complete, reviewed manuscript by **Sunday, October 25**; October 26–30 are reserved for corrections and submission checks.
- Six US-letter IEEE-format pages **including references**; up to two paid extra pages are allowed. Draft for six pages. The call encourages a supplementary video (MP4, at most 3 minutes and 20 MB).
- Double-anonymous review: anonymize manuscript, figures, videos, file metadata, and repository links. Follow the call's AI-content disclosure policy in the submitted version as applicable, while preserving review anonymity.

## Research question and claim

**Question:** Can a compact, physically structured two-DOF dynamic model predict the measured bending of this three-PMA module under time-varying measured chamber pressures, and improve tracking over simpler pressure-to-shape control?

**Candidate contribution, conditional on results:** A calibrated pressure-actuation interface and identified stiffness/damping for the existing single-section dynamics, with held-out open-loop prediction and a controlled tracking experiment on the physical module. Report the tested motion/pressure envelope. Do not claim the first recursive model, the first pressure-coordinate model, or a complete valve-flow model.

**Naming contract:** `u ∈ [0,10]^3 V` is the regulator command; `p_ch ∈ R^3` is **measured physical chamber gauge pressure**; `z ∈ R^2` is independently measured mechanical bending/shape; `q_eq ∈ R^2` is the existing pressure-equivalent *shape coordinate*; `Q ∈ R^2` is generalized effort. The existing `arm_dynamics_pressure_2DoF_MSF.m` uses `tau=Q`, not `u` or `p_ch`. The zero-sum relation used in the shape chart is not a physical constraint on the three nonnegative chamber pressures.

**Model bridge to establish:**

```text
voltage u  -->  measured/regulator pressure p_ch  -->  Q = B(z,p_ch; beta)
                                                     -->  M(q_eq) q_eq_dd + C q_eq_d
                                                          + G + F_el + D q_eq_d = Q
                                                     -->  measured shape z, tip pose
```

Use measured `p_ch(t)` as the mechanical-model input during identification and validation. Model the regulator's command-to-pressure response separately only as needed for the closed-loop controller. The pressure-to-effort mapping `B` must be derived or calibrated and its feasible pressure allocation stated. Avoid fitting an unconstrained scale into both the actuation gain and stiffness without an identifiability check. The current definition of `K` and its units must be reconciled before reporting it as physical bending rigidity; if that is not possible by the paper freeze, report identified **effective** parameters and their coordinate convention.

## Six-page manuscript layout

| Space | Section | Required content and evidence |
| --- | --- | --- |
| Page 1 | Title, abstract, introduction | Robot photo/diagram; problem; precise contribution; one-sentence numerical result placeholders to fill only after experiments; related-work distinction. |
| Page 2 | Module and model | Three-PMA hardware and two bending coordinates; measurement setup; equation of motion; `u`, `p_ch`, `q_eq`, `Q` distinction; assumptions and operating range. |
| Page 3 | Pressure coupling and identification | Actuation map, regulator response if used, excitation design, parameter vector, bounded least-squares objective, training/validation split, identifiability statement. |
| Page 4 | Dynamic-model validation | Held-out measured-vs-predicted trajectories in both bending directions; static-map and simple dynamic baselines; pressure traces; error table and repeat variability. |
| Page 5 | Control demonstration | Model-based feedforward plus measured-shape feedback, feasible pressure allocation, saturation handling; compare with direct inverse pressure map and tuned feedback-only control on the same unseen references. |
| Page 6 | Discussion, limitations, conclusion, references | State where the model succeeds/fails, computational cost, sensitivity to pressure/shape measurement and unmodeled hysteresis; concise conclusion and references. Reserve enough space for references from the first draft. |

**Provisional title:** “Identification and Experimental Validation of a Pressure-Driven Dynamic Model for a Hybrid Soft Continuum Module.” Revise after results; include “control” in the title only if the hardware control comparison is convincing.

**Figure/table budget:** (1) module and sensing schematic; (2) model/input notation and measured regulator response; (3) held-out two-axis trajectory overlay with three measured chamber-pressure traces; (4) control tracking overlay; (5) compact parameter/error/timing table. Combine panels to stay within six pages. Supplementary video should show real hardware and synchronized overlays if available; an animation alone is not experimental validation.

## Minimum executable experiment set

1. **Instrumentation check:** Synchronize timestamps for three voltages, three chamber pressures, and independently measured shape/tip pose. Record sampling rates, sensor calibration, pressure limits, and a zero-pressure baseline. Save raw data and metadata before processing. If chamber pressure or shape is unavailable, do not label a voltage-to-shape fit as mechanical parameter identification.
2. **Regulator characterization:** Use modest steps or chirps to compare commanded and measured pressure, including both inflation and venting. Choose the shape-excitation frequency range after checking the measured pressure bandwidth; do not select a chirp range solely from command voltage.
3. **Identification:** Acquire slow/static holds at several pressure combinations for equilibrium behavior. Excite two independent differential-pressure directions around a feasible positive bias; test common pressure as a separate condition because it may alter length or apparent stiffness. Fit a compact shared parameter set with `lsqnonlin` or equivalent bounded optimization. Use measured pressure as input and forward-simulated shape as residual. Start with measured module mass/inertia and geometric constants; estimate only parameters the data can support. Repeat at least one trial to assess variability.
4. **Held-out validation:** Freeze parameters; predict new waveform(s), amplitude(s), and initial pose(s). Include at least one shape direction not used as a sole-axis identification trial and a return/venting phase. Report RMSE or normalized RMSE in both bending coordinates and tip position, phase/timing error, and visible failure cases. Compare with a quasi-static pressure-to-shape map and a low-order dynamic baseline using the same training data.
5. **Control:** If the calibrated model and feasible pressure map pass validation, run matched reference trajectories with (a) direct static pressure mapping, (b) feedback-only shape control, and (c) model-based feedforward plus the same feedback. Use actual pressure and shape feedback. Report tracking error, settling/phase, peak pressure, saturation, and repeated-trial variation. Do not claim feedback linearization unless the physical input map, allocation, and constraints support it.

**Data split:** Name training trials before fitting; keep validation and control reference trials untouched until parameter freeze. Keep raw data, processing scripts, fit settings, and plots reproducible. A simulation video from `runDynamicSimulation_pressure_2DoF_MSF.m` is illustrative only.

## Schedule to the Sunday manuscript freeze

| Date (Chicago time) | Deliverable and decision |
| --- | --- |
| Mon Oct 19 | Freeze paper question, coordinate/input notation, experimental channels, six-page outline, and figure list. Check that pressure and shape can be recorded together. |
| Tue Oct 20 | Make hardware data logging reliable; record sensor/clock calibration and regulator response. Draft introduction, hardware, and model sections with placeholder-free methods. |
| Wed Oct 21 | Complete identification trials in two bending directions plus static holds; archive raw data. **Gate 1:** if synchronized pressure and shape data are not usable, narrow the claim immediately. |
| Thu Oct 22 | Fit parameters, check residuals and parameter correlation, freeze training data and parameter set. Draft identification section. |
| Fri Oct 23 | Run held-out dynamic validation and baseline comparisons; generate main validation figure/table. **Gate 2:** if prediction does not beat simple baselines, revise the model or paper claim before control work. |
| Sat Oct 24 | Run matched control comparisons if the model is ready; otherwise devote the day to additional validation and a modeling-only paper. Produce control figure/video if supported. |
| Sun Oct 25 | Complete six-page manuscript, references, abstract with measured numbers, internal technical review, anonymization pass, and PDF/figure checks. Freeze a reviewable version. |
| Oct 26–30 | Address review findings, rerun only experiments needed to resolve concrete claims, check PaperPlaza compliance, upload early. Official deadline: Oct 31, 11:59 PM Pacific. |

## Responsibilities and current unknowns

- **Researcher/hardware team:** Physical excitation and safety limits, sensor calibration, synchronized pressure/shape recordings, controller deployment, and experimental sign-off.
- **Codex/code and writing support:** Audit the single-section equations and units; implement reproducible import, identification, simulation, comparison, and plotting scripts once the data format is known; draft and revise the paper from verified results; maintain `PROJECT_LOG.md` with decisions and measured evidence.
- **Need to confirm now:** Available shape sensor or motion capture, actual per-PMA chamber-pressure sensing, logging rate/timestamp synchronization, pressure range, pressure-regulator interface, experiment time, and who can run the physical controller. These determine whether the control result is achievable by October 24.

## Paper integrity gates

- Report estimated **effective** stiffness/damping if physical stiffness and actuator gain cannot be separated. Do not recycle simulation parameters as identified measurements.
- State the pressure-coordinate chart and its dependence on `K`; derive outputs from independent measured shape. Resolve the current `K` dimensional inconsistency before claiming a physical modulus or rigidity.
- State limitations from hysteresis, common pressure, extension, and unmodeled regulator dynamics. Keep all numerical claims traceable to frozen scripts and raw trials.
- If the weekend evidence supports only identification and prediction, submit a modeling/validation paper only if that contribution is sufficiently distinct from prior work; otherwise use the February 2027 non-proceedings extended abstract or a later venue.
