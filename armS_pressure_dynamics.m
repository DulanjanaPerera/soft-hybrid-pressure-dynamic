function dX = armS_pressure_dynamics(~,X,params,Q)
% Interpreted three-section pressure dynamics, X=[q;dq] (12x1).
% M*qdd+C*dq+G+Fel+D*dq=Q. Q is a work-conjugate generalized force
% (J/Pa), not a pressure command or pneumatic regulator model.
assert(isequal(size(X),[12,1]) && isequal(size(Q),[6,1]));
assert(isequal(size(params.D),[6,6]));
assert(all(isfinite(X)) && all(isfinite(Q)) && ...
    norm(params.D-params.D.','fro')<1e-12 && ...
    all(eig(params.D)>=0));
q = X(1:6); dq = X(7:12);
[M,C,G] = armS_pressure_core(q,dq,params);
Fel = pressure_elastic_force(q,params);
ddq = M\(Q-C*dq-G-Fel-params.D*dq);
dX = [dq;ddq];
end
