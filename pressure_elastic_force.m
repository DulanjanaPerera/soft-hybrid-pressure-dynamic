function [Fel,U,H] = pressure_elastic_force(q,params)
% Pressure-coordinate elastic force and potential for three sections.
% q = [p12;p13;p22;p23;p32;p33] (Pa). Fel = dU/dq has units J/Pa.
% This is the elastic part of G_and_K_pressure_2DoF_MSF at g=0, applied
% independently to each section. k(1,n) is not explicit in this formula;
% it can enter through the chosen section bending stiffness K(n).

assert(isequal(size(q),[6,1]));
assert(isequal(size(params.k),[3,3]));
assert(isequal(size(params.r),[3,1]) && ...
    isequal(size(params.K),[3,1]) && isequal(size(params.A),[3,1]));
assert(all(isfinite(q)) && all(params.k(:)>0) && ...
    all(params.r>0) && all(params.K>0) && all(params.A>0));

H = zeros(6);
for n = 1:3
    rows = 2*n-1:2*n;
    s = params.A(n)^2*params.r(n)^2/params.K(n);
    spring = 9*params.r(n)^2/(4*params.K(n));
    H(rows,rows) = s*[3+spring*params.k(2,n), 3/2; ...
                       3/2, 3+spring*params.k(3,n)];
end
Fel = H*q;
U = 0.5*q.'*Fel;
end
