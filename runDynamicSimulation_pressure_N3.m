% Run and record the interpreted three-section pressure-coordinate model.
% Edit the parameter, initial-state, input, and video blocks below.
% State order: X = [q;dq], q = [p12;p13;p22;p23;p32;p33] in Pa.
% Each column of params.L describes one section:
% [base sensor offset; flexible backbone length; tip sensor offset] (m).
% Sensor offsets do not add backbone length or mass to this model.

%% Section parameters (columns 1, 2, and 3)
params.L = [0.001, 0.001, 0.001; ... % base sensor offsets (m)
            0.178, 0.178, 0.178; ... % physical section lengths (m)
            0.001, 0.001, 0.001];    % tip sensor offsets (m)
params.r = [0.013;0.013;0.013];     % actuator radial offsets (m)
params.A = pi*(params.r/2).^2;       % effective actuator areas (m^2)
params.mi = [0.10;0.1;0.10];        % total mass of each section (kg)

%% Base pose and gravity
% Rotate the whole arm by editing these angles (degrees). Rotations are
% applied in yaw-pitch-roll order: Rbase = Rz(yaw)*Ry(pitch)*Rx(roll).
% A 180 degree roll or pitch reverses the backbone's initial direction.
baseRollDeg = 0;
basePitchDeg = 0;
baseYawDeg = 0;
params.basePosition = [0;0;0];      % base location in world coordinates (m)
Rx = [1,0,0;0,cosd(baseRollDeg),-sind(baseRollDeg); ...
    0,sind(baseRollDeg),cosd(baseRollDeg)];
Ry = [cosd(basePitchDeg),0,sind(basePitchDeg);0,1,0; ...
    -sind(basePitchDeg),0,cosd(basePitchDeg)];
Rz = [cosd(baseYawDeg),-sind(baseYawDeg),0; ...
    sind(baseYawDeg),cosd(baseYawDeg),0;0,0,1];
params.baseRotation = Rz*Ry*Rx;
params.gWorld = [0;0;9.81];         % gravitational acceleration in world frame (m/s^2)
params.g = params.baseRotation.'*params.gWorld; % same gravity in base frame

% PMA stiffness: rows are local actuators 1, 2, 3 (N/m).
% The current elastic law uses rows 2 and 3 directly. Row 1 enters the
% example K relation below; adjust K directly if measured values exist.
params.k = 3200*ones(3,3);
params.K = (3/20)*params.k(1,:).'.*params.L(2,:).'.*params.r.^2;
% To supply measured bending stiffnesses instead, replace params.K with
% a positive 3x1 vector. The present model's K units need review.

% Viscous Rayleigh damping in pressure coordinates (6x6, symmetric PSD).
params.D = 8e-10*eye(6);

%% Initial state
q0 = [10000;0; ...  % section 1: p12, p13 (Pa)
       200;100; ... % section 2: p22, p23 (Pa)
       100;0];       % section 3: p32, p33 (Pa)
dq0 = zeros(6,1);     % pressure-coordinate velocities (Pa/s)
X0 = [q0;dq0];

%% Applied generalized force
% This is work-conjugate to q (units J/Pa), NOT a pressure command.
% A constant 6x1 vector or a function @(time,state) returning 6x1 is
% accepted. For example: inputForce = @(time,state) [1e-5*sin(2*pi*time);zeros(5,1)];
inputForce = zeros(6,1);

%% Integration settings
tFinal = 10;                            % s
nSamples = 2001;                         % reported samples, including ends
relativeTolerance = 1e-8;
absoluteTolerance = [1e-5*ones(6,1); ... % Pa
                     1e-4*ones(6,1)];   % Pa/s

%% Animation and video settings
makeAnimation = true;
makePressurePlots = true;             % pressure traces in Figure 2
videoOptions.saveVideo = true;
videoOptions.fileName = fullfile('results','armS_pressure_run.mp4');
videoOptions.frameStride = 8;
videoOptions.fps = 25;
videoOptions.nXi = 31;
videoOptions.nCircle = 16;
videoOptions.bodyRadius = 0.015;       % m
videoOptions.visible = 'on';
videoOptions.closeFigure = false;

%% Run
assert(isequal(size(q0),[6,1]) && isequal(size(dq0),[6,1]));
assert(isscalar(tFinal) && tFinal>0 && nSamples>=2);
assert(isequal(size(absoluteTolerance),[12,1]) && ...
    all(absoluteTolerance>0));
if isa(inputForce,'function_handle')
    Q0 = inputForce(0,X0);
    assert(isequal(size(Q0),[6,1]) && all(isfinite(Q0)));
    rhs = @(time,state) armS_pressure_dynamics( ...
        time,state,params,inputForce(time,state));
else
    assert(isequal(size(inputForce),[6,1]) && all(isfinite(inputForce)));
    rhs = @(time,state) armS_pressure_dynamics( ...
        time,state,params,inputForce);
end
sampleTimes = linspace(0,tFinal,nSamples);
solverOptions = odeset('RelTol',relativeTolerance, ...
    'AbsTol',absoluteTolerance);
tic
[t,X] = ode15s(rhs,sampleTimes,X0,solverOptions);
simulationSeconds = toc;
fprintf('Three-section simulation: %d samples over %.3f s; elapsed %.2f s.\n', ...
    numel(t),t(end),simulationSeconds);

if makePressurePlots
    plot_armS_pressure_history(t,X);
end

videoFile = '';
if makeAnimation
    videoFile = animate_armS_pressure(t,X,params,videoOptions);
    if videoOptions.saveVideo
        fprintf('Animation recorded: %s\n',videoFile);
    end
end
