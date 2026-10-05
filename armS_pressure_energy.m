function [E,T,U,V] = armS_pressure_energy(q,dq,params)
% Mechanical energy: translational kinetic + elastic + gravitational.
% The backbone mass model uses L(2,n); L(1,n), L(3,n) are sensor offsets.
assert(isequal(size(q),[6,1]) && isequal(size(dq),[6,1]));
M = armS_pressure_core(q,dq,params);
T = 0.5*dq.'*M*dq;
[~,U] = pressure_elastic_force(q,params);
V = 0;
Pbase = zeros(3,1); Rbase = eye(3);
for n = 1:3
    p = [0,q(2*n-1:2*n).'];
    Lbody = [0;params.L(2,n);0];
    [mu] = integratedPosition_nume(p,Lbody,params.r(n), ...
        params.K(n),params.A(n));
    V = V-params.mi(n)*params.g.'*(Pbase+Rbase*mu);
    [~,Rtip,ptip] = HTM_nume_mex(p,1,Lbody,params.r(n), ...
        params.K(n),params.A(n));
    Pbase = Pbase+Rbase*ptip;
    Rbase = Rbase*Rtip;
end
E = T+U+V;
end
