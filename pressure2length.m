function l = pressure2length(p, r, K, A)
% Compute the lengths from the pressure. Here only length-2 and length-3 is
% computed. The length-1 is derived from following
% 
%   l1 = -(l2 + l3)
% 
% Inputs:
%   p  : Pressures of 2nd and 3rd PMAs [1x3] (Pa)
%   r  : radial offset of the PMA from the centerline [constant] (m)
%   K  : Bending stiffness of the module [consatant] (N/rad)
%   A  : PMA effective area [consatant] (m^2)
% 
% Output:
%   l  : l vector (only l2 and l3) [2x1] (m)


arguments (Input)
    p (3,1) double 
    r (1,1) double {mustBePositive} = 0.013
    K (1,1) double {mustBePositive} = 0.031
    A (1,1) double {mustBePositive} = pi*(0.013/2)^2
end

arguments (Output)
    l (2,1) double
end


l = [0.8660254040e0 * r ^ 2 * A * p(2) / K * (p(2) ^ 2 + p(2) * p(3) + p(3) ^ 2) ^ (-0.1e1 / 0.2e1) * sqrt(0.3e1 * p(2) ^ 2 + 0.3e1 * p(2) * p(3) + 0.3e1 * p(3) ^ 2) 0.8660254040e0 * r ^ 2 * A * p(3) / K * (p(2) ^ 2 + p(2) * p(3) + p(3) ^ 2) ^ (-0.1e1 / 0.2e1) * sqrt(0.3e1 * p(2) ^ 2 + 0.3e1 * p(2) * p(3) + 0.3e1 * p(3) ^ 2)];

l = reshape(l, [2,1]);

end