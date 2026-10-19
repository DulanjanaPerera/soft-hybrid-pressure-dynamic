# RoboSoft 2027 single-module logging specification

Status: acquisition contract, October 19, 2026. Available measured channels are the **three commanded regulator pressure setpoints**, tip XYZ, and `theta,phi`. Actual PMA chamber pressures are **not** measured.

## Raw trial files

Create one directory per trial, for example `trial_001/`. Preserve the unmodified acquisition files. Use a shared monotonic time base, in seconds, if commands and pose are produced by the same computer. Do not replace raw data with interpolated or filtered traces.

Preferred format for separate command and pose acquisition loops:

`commands.csv`

```csv
t_s,p_cmd1_Pa,p_cmd2_Pa,p_cmd3_Pa
0.000,0,0,0
```

`pose.csv`

```csv
t_s,tip_x_m,tip_y_m,tip_z_m,theta_rad,phi_rad
0.003,0.0,0.0,0.278,0.0,0.0
```

Each `commands.csv` row records the setpoints **actually sent** and the time they were sent, including initial and final commands. Each `pose.csv` row records the most recent pose measurement and its acquisition time. The example numbers above are format examples, **not experimental data**. Keep time increasing in each file. Do not silently duplicate a stale pose value when a measurement is missing; use a missing-value marker and record it in metadata.

If all channels are sampled in one synchronized loop, one file `samples.csv` is acceptable with `t_s` and all eight signal columns above. If commands and pose have independent clocks, keep their original device timestamps and record the clock-sync method and observed offset; do not pretend their timestamps share an origin.

## Coordinate and units metadata

Add `trial.json` or a short text file with:

- Trial ID, date/time, experimenter, module identifier, firmware/software version.
- Waveform type (static hold, step, chirp, combined-axis validation, or control), which channels were excited, commanded bias/amplitude/frequency range, and intended training/validation/control split **assigned before fitting**.
- Commanded-pressure unit and whether commands are gauge pressure setpoints. Convert bar or kPa to Pa in the analysis copy while retaining original raw units in the acquisition file if needed. Record regulator command limits and any software saturation.
- Tip XYZ reference frame, origin, axis directions, and whether the reported pose is measured directly or computed from the same markers/sensor as `theta,phi`.
- Exact definitions of `theta` and `phi`, angle unit, sign, range, and wrap convention. The repository's `pressure2config.m` uses its first output as a bending-plane angle and second as bend magnitude; check the experimental logger against this before comparing trajectories. A bending-plane angle becomes ill-defined near the straight pose.
- Pose acquisition rate, command update rate, known sensor latency or buffering, missing-sample conventions, and any filtering performed before logging.
- Physical setup: base orientation, module length, attached payload if any, supply setting, and any change to PMA plumbing or fittings.

## One-minute commissioning trial

1. With the module in a safe state, record a short zero-command baseline and verify the tip pose is finite and stable.
2. Command one modest, safe step in PMA 1; return to baseline. Repeat for PMAs 2 and 3. Verify timestamps, command order, sign/direction of tip motion, return behavior, and that no channel is swapped.
3. Inspect elapsed duration, command/pose sampling gaps, dropped pose samples, and possible time offset. Plot commanded pressure and tip XYZ/angles versus time; do not infer unmeasured chamber pressure from these plots.
4. Freeze the raw file and note any issue before collecting identification chirps. Choose chirp frequency range only after viewing this command-to-pose response and the sensor sampling rate.

## Analysis contract

Use `theta,phi` as primary dynamic-fit outputs. Predict tip XYZ through the model geometry and report its error as an additional output check. If XYZ and angles are computed from the same sensor data, they are correlated outputs and should not be described as independent measurements. During preprocessing, align commands and pose using documented timestamps, retain raw copies, and record any interpolation/filtering and angle unwrapping. Use a training/validation split at the **trial** level, never random time samples from the same trajectory.

The fitted model is a command-to-motion model. Without actual chamber-pressure or force sensing, fitted stiffness, damping, and actuator gain remain **effective** parameters under the regulator/PMA response and chosen coordinate chart. A release/free-decay trial can constrain the response but cannot prove immediate venting of the chambers.

## MATLAB/Simulink export

The helper [save_robosoft_trial.m](save_robosoft_trial.m) accepts timestamped MATLAB tables and writes the two CSVs and `trial.json`. It checks required columns, strictly increasing timestamps, pressure units, and existing-file protection. It does **not** acquire hardware data or synchronize separate device clocks.

For a Simulink model, log a 3-element vector of the **actual command values sent** and a 5-element vector `[tip_x;tip_y;tip_z;theta;phi]` with `To Workspace` blocks in `Timeseries` format. If both have the same trustworthy acquisition clock and the exported variables are `pCmdTS` and `poseTS`, a conversion looks like:

```matlab
tc = pCmdTS.Time(:);
pc = reshape(pCmdTS.Data, numel(tc), 3);
commands = table(tc,pc(:,1),pc(:,2),pc(:,3), ...
    'VariableNames',{'t_s','p_cmd1_Pa','p_cmd2_Pa','p_cmd3_Pa'});

tp = poseTS.Time(:);
yp = reshape(poseTS.Data, numel(tp), 5);
pose = table(tp,yp(:,1),yp(:,2),yp(:,3),yp(:,4),yp(:,5), ...
    'VariableNames',{'t_s','tip_x_m','tip_y_m','tip_z_m', ...
                     'theta_rad','phi_rad'});

metadata = struct('trial_id','trial_001', ...
    'trial_role','commissioning', ...
    'command_unit','Pa', ...
    'pose_frame','describe base frame here', ...
    'theta_definition','describe angle and sign here', ...
    'phi_definition','describe angle and sign here');
save_robosoft_trial("trial_001",commands,pose,metadata);
```

Replace placeholder definitions with the measured conventions. If the command uses bar or kPa, convert the *analysis/export copy* to Pa before constructing `commands`, while preserving the original raw log. `Timeseries.Time` may be simulation time rather than the real pose-acquisition instant in hardware-in-the-loop or asynchronous sensing; in that case, log hardware timestamps instead. MATLAB can create the same two tables directly from a live acquisition loop, with one timestamp taken when the command is dispatched and one when each pose is acquired.
