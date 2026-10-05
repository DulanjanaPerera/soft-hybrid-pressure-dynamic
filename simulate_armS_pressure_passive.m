function [t,X,params,diagnostics] = simulate_armS_pressure_passive()
% Interpreted passive release of a three-section pressure-coordinate arm.
% The initial pressures are generalized coordinates, not actuator commands.
params.L = [0.010,0.008,0.012; ...
            0.278,0.260,0.290; ...
            0.015,0.012,0.006];
params.r = [0.013;0.0125;0.014];
params.k = 3200*ones(3,3);
params.K = (3/20)*params.k(1,:).'.*params.L(2,:).'.*params.r.^2;
params.A = pi*(params.r/2).^2;
params.mi = [0.10;0.08;0.12];
params.g = [0;0;-9.81];
params.D = 8e-11*eye(6);
Q = zeros(6,1);
X0 = [[10000;0;2000;1000;1000;0];zeros(6,1)];
sampleTimes = linspace(0,2,401);
options = odeset('RelTol',1e-8,'AbsTol',[1e-5*ones(6,1); ...
    1e-4*ones(6,1)]);
[t,X] = ode15s(@(time,state) armS_pressure_dynamics( ...
    time,state,params,Q),sampleTimes,X0,options);

diagnostics.energy = zeros(numel(t),1);
diagnostics.kinetic = zeros(numel(t),1);
diagnostics.elastic = zeros(numel(t),1);
diagnostics.gravity = zeros(numel(t),1);
diagnostics.dissipationRate = zeros(numel(t),1);
for j = 1:numel(t)
    q = X(j,1:6).'; dq = X(j,7:12).';
    [diagnostics.energy(j),diagnostics.kinetic(j), ...
        diagnostics.elastic(j),diagnostics.gravity(j)] = ...
        armS_pressure_energy(q,dq,params);
    diagnostics.dissipationRate(j) = dq.'*params.D*dq;
end
diagnostics.energyBalanceResidual = diagnostics.energy(end) ...
    -diagnostics.energy(1)+trapz(t,diagnostics.dissipationRate);
diagnostics.maxPressure = max(abs(X(:,1:6)),[],1);
end
