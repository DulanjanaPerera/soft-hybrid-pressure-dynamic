function P = f20251121_6_backbonePos_pressure_xi_2DoF_MSF(p, L, r, xi, K, A)
% The backbone point at xi
% Inputs:
%   p  : pressures of the PMAs [1x3] (Pa) (however, p(1) = 0)
%   L  : neutral length of the backbone [3x1] (m)
%   r  : radial offset of the PMA from the centerline [constant] (m)
%   xi : selection parameter on the backbone [constant {0-1}]
%   K  : Module bending stiffness [constant] (N/rad)
%   A  : Effective aream of the PMA [constant] (m^2)
% 
% Output:
%   P  : position vector [3x1]

arguments (Input)
    p (3,1) double 
    L (3,1) double {mustBeGreaterThanOrEqual(L,0)} = [0, 0.278, 0]
    r (1,1) double {mustBePositive} = 0.013
    xi (1,1) double {mustBeBetween(xi, 0, 1, "closed")} = 1
    K (1,1) double {mustBePositive} = 0.031
    A (1,1) double {mustBePositive} = pi*(0.013/2)^2
end

arguments (Output)
    P (3,1) double
end


P = [-0.9e1 / 0.8960e4 * (p(2) + p(3)) * xi ^ 2 * A * L(2) * r * (xi ^ 6 * A ^ 6 * r ^ 6 * (p(2) ^ 2 + p(2) * p(3) + p(3) ^ 2) ^ 3 - 0.56e2 / 0.3e1 * xi ^ 4 * A ^ 4 * r ^ 4 * (p(2) ^ 2 + p(2) * p(3) + p(3) ^ 2) ^ 2 * K ^ 2 + 0.560e3 / 0.3e1 * xi ^ 2 * A ^ 2 * r ^ 2 * (p(2) ^ 2 + p(2) * p(3) + p(3) ^ 2) * K ^ 4 - 0.2240e4 / 0.3e1 * K ^ 6) / K ^ 7 0.3e1 / 0.8960e4 * (p(2) - p(3)) * xi ^ 2 * A * L(2) * r * sqrt(0.3e1) * (xi ^ 6 * A ^ 6 * r ^ 6 * (p(2) ^ 2 + p(2) * p(3) + p(3) ^ 2) ^ 3 - 0.56e2 / 0.3e1 * xi ^ 4 * A ^ 4 * r ^ 4 * (p(2) ^ 2 + p(2) * p(3) + p(3) ^ 2) ^ 2 * K ^ 2 + 0.560e3 / 0.3e1 * xi ^ 2 * A ^ 2 * r ^ 2 * (p(2) ^ 2 + p(2) * p(3) + p(3) ^ 2) * K ^ 4 - 0.2240e4 / 0.3e1 * K ^ 6) / K ^ 7 (xi ^ 8 * A ^ 8 * r ^ 8 * (p(2) ^ 2 + p(2) * p(3) + p(3) ^ 2) ^ 4 - 0.24e2 * xi ^ 6 * A ^ 6 * r ^ 6 * (p(2) ^ 2 + p(2) * p(3) + p(3) ^ 2) ^ 3 * K ^ 2 + 0.336e3 * xi ^ 4 * A ^ 4 * r ^ 4 * (p(2) ^ 2 + p(2) * p(3) + p(3) ^ 2) ^ 2 * K ^ 4 - 0.2240e4 * xi ^ 2 * A ^ 2 * r ^ 2 * (p(2) ^ 2 + p(2) * p(3) + p(3) ^ 2) * K ^ 6 + 0.4480e4 * K ^ 8) * xi * L(2) / K ^ 8 / 0.4480e4];

P = reshape(P, [3,1]);

end
