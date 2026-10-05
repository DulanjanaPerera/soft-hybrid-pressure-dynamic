function validate_pressure_local_interfaces()
% Check the fixed input and output contracts of the pressure-space exports.
% Local coordinates are [0,p2,p3] in Pa. Sensor offsets are not mass.

p = [0, 2000, 1000];
L = [0.01; 0.278; 0.015];
Lbody = [0; L(2); 0];
r = 0.013;
K = 0.0226;
A = pi*(0.013/2)^2;
xi = 0.4;

[T,R,P] = HTM_nume_mex(p,xi,L,r,K,A);
[PosJ,RotJ,PosJJ,RotJJ] = LocalJacob_nume_mex(p,xi,L,r,K,A);
assert(isequal(size(T),[4,4]) && isequal(size(R),[3,3]) ...
    && isequal(size(P),[3,1]));
assert(isequal(size(PosJ),[3,2]) && isequal(size(RotJ),[3,6]) ...
    && isequal(size(PosJJ),[6,2]) && isequal(size(RotJJ),[6,6]));
assert(all(isfinite(T(:))) && all(isfinite(PosJ(:))) ...
    && all(isfinite(RotJ(:))) && all(isfinite(PosJJ(:))) ...
    && all(isfinite(RotJJ(:))));
assert(norm(T(1:3,1:3)-R,'fro') < 1e-12);
assert(norm(T(1:3,4)-P) < 1e-12);

straight = [0,0,0];
[~,~,Pbase] = HTM_nume_mex(straight,0,Lbody,r,K,A);
[~,~,Ptip] = HTM_nume_mex(straight,1,Lbody,r,K,A);
assert(norm(Pbase) < 1e-12);
assert(norm(Ptip-[0;0;L(2)]) < 1e-10);

% Mass integrals must use the flexible section length without sensor offsets.
[mu,mu_q,mu_qq] = integratedPosition_nume(p,Lbody,r,K,A);
[S,S_q] = integratedPositionProduct_nume(p,Lbody,r,K,A);
[F,F_q] = integratedPositionDerivativeProduct_nume(p,Lbody,r,K,A);
[E,E_q] = integratedJacobianProduct_nume(p,Lbody,r,K,A);
values = {mu,mu_q,mu_qq,S,S_q,F,F_q,E,E_q};
sizes = {[3,1],[3,2],[3,2,2],[3,3],[3,3,2], ...
    [3,3,2],[3,3,2,2],[2,2],[2,2,2]};
for i = 1:numel(values)
    assert(isequal(size(values{i}),sizes{i}), ...
        'Unexpected output size at item %d.',i);
    assert(all(isfinite(values{i}(:))), ...
        'Nonfinite output at item %d.',i);
end
[muStraight,~,~] = integratedPosition_nume(straight,Lbody,r,K,A);
assert(norm(muStraight-[0;0;L(2)/2]) < 1e-10);

% A is an explicit runtime input, not a value from the caller workspace.
[~,~,PsmallA] = HTM_nume_mex(p,xi,L,r,K,A/2);
assert(norm(PsmallA-P) > 1e-8);
fprintf('Pressure local function interfaces passed.\n');
end
