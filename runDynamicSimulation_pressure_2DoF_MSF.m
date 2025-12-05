% This is the MATLAB script that run the dynamic simulation of a soft
% continuum arm. The continuum arm is hybrid design which contains both
% soft and rigid elements to achieve deformation and the structural
% rigidity. Please refer following research articles for more informations.
% 
%  Research Articles:
%   [1]   Arachchige, Dimuthu DK, Yue Chen, Ian D. Walker, and Isuru S. Godage. 
%       "A novel variable stiffness soft robotic gripper." In 2021 IEEE 17th 
%       International Conference on Automation Science and Engineering (CASE), 
%       pp. 2222-2227. IEEE, 2021.
%   [2]   Arachchige, Dimuthu DK, and Isuru S. Godage. "Hybrid soft robots 
%       incorporating soft and stiff elements." In 2022 IEEE 5th international 
%       conference on soft robotics (RoboSoft), pp. 267-272. IEEE, 2022.
% 
% The dynamic model is derived as presented in the following article
%   [3]   Godage, Isuru S., David T. Branson, Emanuele Guglielmino, 
%       Gustavo A. Medrano-Cerda, and Darwin G. Caldwell. "Dynamics for biomimetic 
%       continuum arms: A modal approach." In 2011 IEEE International Conference on 
%       Robotics and Biomimetics, pp. 104-109. IEEE, 2011.
% 
% However, the dynamic model is derived for actuator-space variables and the input to
% the system is control signal of the pressure regulator. Since, backbone
% contraint, only two actuator degress are required for the task-space
% control. There for out of 3 PMAs pressures, ob=nly two is required in the
% model.
%
%   p1 = -(p2 + p3); This is the relationship of the actuator-space
% 
% The equation of motion is in the form of:
%       M(l)*p_ddot + C(p, p_dot)*p_dot + G(p) = V
%           
%   p       - 2 pressures of PMAs (lenght of each actuator) (m) [2x1]
%   p_dot   - pressure speed (m/s) [2x1]
%   p_ddot  - pressure acceleration (m/s^2) [2x1]
%   signal  - Pressures (bar) [2x1]
% 
% NOTE: The modal form is used to remove the singularities [3].
% 
% The parameters of the continuum arm is defined with ISO units. The most
% of the parameters are for pneumatic muscle actuators (PMAs). The other
% paramerters are for the in-extensible backbone

clear

params.m = 0.1; % mass of the module in kg
params.L = [0.0001, 0.278, 0.0001]; % module backbone length in m
params.A = pi*(0.013/2)^2; % area of the pneumatic fitting in m^2
params.r = 0.013; % radial offset from the backbone in m
params.k = 3.2e3*[1, 1, 1]; % stiffness of each PMA
params.K = 1e-1*(3.0/2)*params.k(1) * params.L(2)*params.r^2; % bending stiffness of the complete module
params.kmax = 1e6 * ones(1,3);
params.lmin = -0.02 * ones(1,3);
params.lmax = 0.02 * ones(1,3);
params.g = 9.81;
params.d = 5e-11*diag([1 1 1]); % original 5e-12
% params.tau = [params.A*3e5; params.A*3e5];
params.tau = [0.0; 0.0];

xi = linspace(0,1,20);

% Initial state
p0  = [10000; 0.0];     % initial joint coords
dp0 = [0.0; 0.0];
X0  = [p0; dp0];

% Integrate
tspan = [0 10];
tic
[t, X] = ode15s(@(t,X) arm_dynamics_pressure_2DoF_MSF(t,X,params), tspan, X0);
toc

% drawing the arm
tip_pos = zeros(length(t), length(xi), 3);
l = zeros(1,3);

figure(1); clf
ax = axes;
h_backbone = plot3(NaN,NaN,NaN,'LineWidth',2); hold on
h_tip      = plot3(NaN,NaN,NaN,'o','MarkerFaceColor',[.8 .2 .2],'MarkerEdgeColor','none');
xlabel('X'); ylabel('Y'); zlabel('Z');
rotate3d 'on'; 
ht = title(ax, 'Continuum arm backbone (current frame) and tip');
grid on; axis equal

for k = 1:length(t)
    p = zeros(3,1);
    p(2:3) = X(k,1:2)';  % update DOFs from your state

    % fill positions along the backbone for all xi
    for j = 1:length(xi)
        pos = backbonePos_pressure_xi_2DoF_MSF(p, params.L, params.r, xi(j), params.K, params.A); % [1x3] or [3x1]
        tip_pos(k,j,1) = pos(1);
        tip_pos(k,j,2) = pos(2);
        tip_pos(k,j,3) = pos(3);
    end

    % extract this frame's curve
    Xs = squeeze(tip_pos(k,:,1));
    Ys = squeeze(tip_pos(k,:,2));
    Zs = squeeze(tip_pos(k,:,3));

    % update plot
    set(h_backbone,'XData',Xs,'YData',Ys,'ZData',Zs);
    set(h_tip,'XData',Xs(end),'YData',Ys(end),'ZData',Zs(end));  % tip = xi(end)
    xlim([-0.3 0.3]); ylim([-0.3 0.3]); zlim([-0.1 0.3]); 
    set(ht, 'String', sprintf('Continuum arm — frame k = %d / %d   t = %.3f s', ...
                          k, numel(t), t(k)));
    drawnow
    pause(0.01);
end

