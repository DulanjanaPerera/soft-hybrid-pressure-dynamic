function [p_out, error_c] = config2pressure(theta, phi, r, K, A)
% The optimization to obtain minimum energy solutions of pressure and
% positive pressure. The function pressure2config.m file compute the
% configuration variable for given pressures.
% 
% Inputs:
%   theta   : desired theta angle [constant] (rad)
%   phi     : desired phi angle [constant] (rad)
%   r       : radial offset of the PMA [constant] (m)
%   K       : bending stffness of the module [constant] (N/rad)
%   A       : effective area of the PMA [constant] (m^2)
% 
% Output:
%   p_est   : estimated positive pressure vector [3x1] (Pa)
% 


    arguments (Input)
        theta (1,1) double {mustBeBetween(theta, -3.1416, 3.1416, "openleft")}
        phi (1,1) double {mustBeBetween(phi, 0, 3.1416, "openleft")}
        r (1,1) double {mustBePositive} = 0.013
        K (1,1) double {mustBePositive} = 0.031
        A (1,1) double {mustBePositive} = pi*(0.013/2)^2
    end
    
    arguments (Output)
        p_out (3,1) double
        error_c (2,1) double
    end

    % initial values for the optimization
    p0 = [0.01; 0.01; 0.001];

    % change the algorithm to sqp
    options = optimoptions('fmincon','Display','off','Algorithm','sqp');

    % defining the cost function and external variables of the function
    fun = @(p) cost(p, theta, phi, r, K, A);

    % Perform the optimization to estimate the pressure
    % The boundary conditions are p = [0, 300000] Pa
    p_est = fmincon(fun, p0, [], [], [], [], [0; 0; 0], [3e5; 3e5; 3e5], [], options);
    
    % Ensure the output pressure vector is positive
    p_est = max(p_est, 0);
    p_out = [-(p_est(2)+p_est(3)); p_est(2); p_est(3)];

    error_c = [theta; phi] - pressure2config(p_out, r, K, A);
end

function c = cost(p, theta, phi, r, K, A)

    % avoiding the singularity
    if phi < 0.0001
        phi = 0.0001;
    end

    config = pressure2config(p, r, K, A);
    c = norm(config - [theta; phi]) + 0.001*norm(p); % adds the regularization of pressure

end