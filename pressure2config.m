function c = pressure2config(p, r, K, A)
% Compute the configuration parameters from the pressure.
% 
%  c = [theta; phi]
% 
% Inputs:
%   p  : Pressures of 2nd and 3rd PMAs [1x3] (Pa)
%   r  : radial offset of the PMA from the centerline [constant] (m)
%   K  : Bending stiffness of the module [consatant] (N/rad)
%   A  : PMA effective area [consatant] (m^2)
% 
% Output:
%   c  : configuration vector (theta, phi) [2x1]


arguments (Input)
    p (3,1) double 
    r (1,1) double {mustBePositive} = 0.013
    K (1,1) double {mustBePositive} = 0.031
    A (1,1) double {mustBePositive} = pi*(0.013/2)^2
end

arguments (Output)
    c (2,1) double
end

c = [atan2(-A * r * sqrt(0.3e1) * (p(2) - p(3)) / 0.2e1, 0.3e1 / 0.2e1 * A * r * (p(2) + p(3))) sqrt(0.3e1) * A * r * sqrt(p(2) ^ 2 + p(2) * p(3) + p(3) ^ 2) / K];

c = reshape(c,[2,1]);  

end