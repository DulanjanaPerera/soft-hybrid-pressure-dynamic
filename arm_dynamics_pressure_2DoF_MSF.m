function dX = arm_dynamics_pressure_2DoF_MSF(t, X, params)
% simulating the arm dynamic with mode shape functions
% Inputs:
%   t : time
%   X : length and length velocity [4x1]
%       [p2, p3, dp2, dp3]
%   params : parameters for matrices (structure)

    

    p = zeros(3,1); % p has to be (3x1) because in the function it is a vector
    p(2:3) = X(1:2);
    % p = min(p, 300000);
    % p = max(p, -300000);

    dp = zeros(3,1); % this is (3x1). because we don't use it in functions
    dp(2:3) = X(3:4);

    m = params.m; % Extracting mass parameters from input
    L = params.L; % Extracting length parameters from input
    r = params.r; % Extracting radius parameters from input
    K = params.K; % bending stiffness
    A = params.A; % get the area
    k = params.k; % Extracting stiffness of PMAs from input
    g = params.g; % Extracting gravitational parameters from input
    d = params.d; % diag([100 100 100]);
    tau = params.tau; % pressure force
    
    dig_val = diag(params.d);
    % Dynamics
    M = M_pressure_2DoF_MSF(p, m, L, r, K, A);
    C = C_pressure_2DoF_MSF(p, m, L, r, K, A, dp);
    D = [dig_val(2)*dp(2); dig_val(3)*dp(3)];
    G = G_and_K_pressure_2DoF_MSF(p, m, L, r, K, A, k, g);


    % temp = tau - C*[dp(2); dp(3)] - D - G;
    % inM = pinv(M);
    ddp = M \ (tau - C*[dp(2); dp(3)] - D - G);

    if any(isnan(ddp), 'all') || any(isinf(ddp), 'all')
        ddp = [0.0; 0.0];
    end
    dX = [dp(2); dp(3); ddp];

end