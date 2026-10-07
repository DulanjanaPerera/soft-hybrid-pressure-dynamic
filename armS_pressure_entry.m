function dX = armS_pressure_entry(t,X,L,r,K,A,mi,g,k,D,Q)
%#codegen
% Fixed-size runtime interface for the three-section pressure dynamics MEX.
% L(:,n) = [sensor offset; backbone length; sensor offset].
params.L = L;
params.r = r;
params.K = K;
params.A = A;
params.mi = mi;
params.g = g;
params.k = k;
params.D = D;
% t is retained for the ODE RHS signature and time-varying caller input.
q = X(1:6);
dq = X(7:12);
[M,C,G] = armS_pressure_core_mex(q,dq,params);
Fel = pressure_elastic_force(q,params);
ddq = M\(Q-C*dq-G-Fel-D*dq);
dX = [dq;ddq];
end
