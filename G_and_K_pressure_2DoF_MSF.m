function G = f20251121_5_G_and_K_pressure_2DoF_MSF(p, m, L, r, K, A, k, g)
% Generalized force matrix and stiffness matrix in jointspace
% Inputs:
%   p  : pressures of the PMAs [1x3] (Pa)
%   m  : mass of the module [constant] (kg)
%   L  : neutral length of the backbone [1x3] (m)
%   r  : radial offset of the PMA from the centerline [constant] (m)
%   K  : Bending stiffness of the module [consatant] (N/rad)
%   A  : PMA effective area [constant] (m^2)
%   k  : stiffness vector of the PMAs [1x3] (N/m)
%   g  : gravity [constant] (m/s^2)
% 
% Output:
%   G  : Gravity matrix [2x1]

arguments (Input)
    p (3,1) double 
    m (1,1) double {mustBePositive}  = 0.1
    L (3,1) double {mustBeGreaterThanOrEqual(L,0)} = [0, 0.278, 0]
    r (1,1) double {mustBePositive} = 0.013
    K (1,1) double {mustBePositive} = 0.031
    A (1,1) double {mustBePositive} = pi*(0.013/2)^2
    k (3,1) double {mustBePositive} = [2.2e3; 2.2e3; 2.2e3]
    g (1,1) double = 9.81
end

arguments (Output)
    G (2,1) double
end

G = [0.1785714286e-3 * A ^ 2 * r ^ 2 * (((0.3499999999e1 * p(2) ^ 6 * p(3) + 0.7499999997e1 * p(2) ^ 5 * p(3) ^ 2 + 0.1000000000e2 * p(2) ^ 4 * p(3) ^ 3 + 0.9500000002e1 * p(2) ^ 3 * p(3) ^ 4 + 0.5999999997e1 * p(2) ^ 2 * p(3) ^ 5 + 0.2499999999e1 * p(2) * p(3) ^ 6 + 0.4999999999e0 * p(3) ^ 7 + 0.1000000000e1 * p(2) ^ 7) * A ^ 6 * r ^ 6 + (-0.2250000000e2 * K ^ 2 * p(2) ^ 5 - 0.5625000001e2 * K ^ 2 * p(2) ^ 4 * p(3) - 0.8999999998e2 * K ^ 2 * p(2) ^ 3 * p(3) ^ 2 - 0.7874999999e2 * K ^ 2 * p(2) ^ 2 * p(3) ^ 3 - 0.4499999999e2 * K ^ 2 * p(2) * p(3) ^ 4 - 0.1125000000e2 * K ^ 2 * p(3) ^ 5) * A ^ 4 * r ^ 4 + (0.2800000000e3 * K ^ 4 * p(2) ^ 3 + 0.4199999999e3 * K ^ 4 * p(2) * p(3) ^ 2 + 0.1400000000e3 * K ^ 4 * p(3) ^ 3 + 0.4199999999e3 * K ^ 4 * p(2) ^ 2 * p(3)) * A ^ 2 * r ^ 2 - 0.1400000000e4 * K ^ 6 * p(2) - 0.6999999999e3 * K ^ 6 * p(3)) * m * g * L(2) + 0.1260000001e5 * p(2) * k(2) * K ^ 6 * r ^ 2 + 0.8399999998e4 * K ^ 7 * p(3) + 0.1680000000e5 * K ^ 7 * p(2)) / K ^ 8; 0.8928571429e-4 * A ^ 2 * r ^ 2 * (((0.1000000000e1 * p(2) ^ 7 + 0.5000000000e1 * p(2) ^ 6 * p(3) + 0.1200000000e2 * p(2) ^ 5 * p(3) ^ 2 + 0.1900000001e2 * p(2) ^ 4 * p(3) ^ 3 + 0.2000000000e2 * p(2) ^ 3 * p(3) ^ 4 + 0.1500000000e2 * p(2) ^ 2 * p(3) ^ 5 + 0.7000000000e1 * p(2) * p(3) ^ 6 + 0.2000000000e1 * p(3) ^ 7) * A ^ 6 * r ^ 6 + (-0.4500000000e2 * K ^ 2 * p(3) ^ 5 - 0.2250000000e2 * K ^ 2 * p(2) ^ 5 - 0.9000000000e2 * K ^ 2 * p(2) ^ 4 * p(3) - 0.1575000000e3 * K ^ 2 * p(2) ^ 3 * p(3) ^ 2 - 0.1800000000e3 * K ^ 2 * p(2) ^ 2 * p(3) ^ 3 - 0.1125000000e3 * K ^ 2 * p(2) * p(3) ^ 4) * A ^ 4 * r ^ 4 + (0.5600000000e3 * K ^ 4 * p(3) ^ 3 + 0.8400000000e3 * K ^ 4 * p(2) ^ 2 * p(3) + 0.8400000000e3 * K ^ 4 * p(2) * p(3) ^ 2 + 0.2800000000e3 * K ^ 4 * p(2) ^ 3) * A ^ 2 * r ^ 2 - 0.2800000000e4 * K ^ 6 * p(3) - 0.1400000000e4 * K ^ 6 * p(2)) * m * g * L(2) + 0.2520000002e5 * p(3) * k(3) * K ^ 6 * r ^ 2 + 0.1680000000e5 * K ^ 7 * p(2) + 0.3360000000e5 * K ^ 7 * p(3)) / K ^ 8;];

G = reshape(G, [2,1]);

end