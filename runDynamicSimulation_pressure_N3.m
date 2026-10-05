% Run and record the interpreted three-section pressure-coordinate model.
% Edit the parameter, initial-state, input, and video blocks below.
% State order: X = [q;dq], q = [p12;p13;p22;p23;p32;p33] in Pa.
% Each column of params.L describes one section:
% [base sensor offset; flexible backbone length; tip sensor offset] (m).
% Sensor offsets do not add backbone length or mass to this model.

%% Section parameters (columns 1, 2, and 3)
params.L = [0.010, 0.008, 0.012; ... % base sensor offsets (m)
            0.278, 0.260, 0.290; ... % physical section lengths (m)
            0.015, 0.012, 0.006];    % tip sensor offsets (m)
params.r = [0.013;0.0125;0.014];     % actuator radial offsets (m)
params.A = pi*(params.r/2).^2;       % effective actuator areas (m^2)
params.mi = [0.10;0.08;0.12];        % total mass of each section (kg)
params.g = [0;0;-9.81];              % gravitational acceleration (m/s^2)

% PMA stiffness: rows are local actuators 1, 2, 3 (N/m).
% The current elastic law uses rows 2 and 3 directly. Row 1 enters the
% example K relation below; adjust K directly if measured values exist.
params.k = 3200*ones(3,3);
params.K = (3/20)*params.k(1,:).'.*params.L(2,:).'.*params.r.^2;
% To supply measured bending stiffnesses instead, replace params.K with
% a positive 3x1 vector. The present model's K units need review.

% Viscous Rayleigh damping in pressure coordinates (6x6, symmetric PSD).
params.D = 8e-11*eye(6);

%% Initial state
q0 = [10000;0; ...  % section 1: p12, p13 (Pa)
       2000;1000; ... % section 2: p22, p23 (Pa)
       1000;0];       % section 3: p32, p33 (Pa)
dq0 = zeros(6,1);     % pressure-coordinate velocities (Pa/s)
X0 = [q0;dq0];

%% Applied generalized force
% This is work-conjugate to q (units J/Pa), NOT a pressure command.
% A constant 6x1 vector or a function @(time,state) returning 6x1 is
% accepted. For example: inputForce = @(time,state) [1e-5*sin(2*pi*time);zeros(5,1)];
inputForce = zeros(6,1);

%% Integration settings
tFinal = 2;                            % s
nSamples = 401;                         % reported samples, including ends
relativeTolerance = 1e-8;
absoluteTolerance = [1e-5*ones(6,1); ... % Pa
                     1e-4*ones(6,1)];   % Pa/s

%% Animation and video settings
makeAnimation = true;
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

videoFile = '';
if makeAnimation
    videoFile = animate_armS_pressure(t,X,params,videoOptions);
    if videoOptions.saveVideo
        fprintf('Animation recorded: %s\n',videoFile);
    end
end
