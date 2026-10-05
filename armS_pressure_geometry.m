function [P,R,bases] = armS_pressure_geometry(q,params,xi)
% Global backbone points/frames for three sections in pressure coordinates.
% P is nXi-by-3-by-3; R is 3-by-3-by-nXi-by-3; bases is 3-by-4.
% Sensor offsets L(1,:) and L(3,:) do not extend the physical backbone.
assert(isequal(size(q),[6,1]) && isvector(xi) && ...
    all(isfinite(xi)) && all(xi>=0) && all(xi<=1));
assert(isequal(size(params.L),[3,3]) && ...
    isequal(size(params.r),[3,1]) && ...
    isequal(size(params.K),[3,1]) && ...
    isequal(size(params.A),[3,1]));
xi = xi(:);
nXi = numel(xi);
P = zeros(nXi,3,3);
R = zeros(3,3,nXi,3);
bases = zeros(3,4);
Pbase = zeros(3,1);
Rbase = eye(3);
for n = 1:3
    p = [0,q(2*n-1:2*n).'];
    Lbody = [0;params.L(2,n);0];
    bases(:,n) = Pbase;
    for j = 1:nXi
        [~,Rlocal,Plocal] = HTM_nume_mex(p,xi(j),Lbody, ...
            params.r(n),params.K(n),params.A(n));
        P(j,:,n) = (Pbase+Rbase*Plocal).';
        R(:,:,j,n) = Rbase*Rlocal;
    end
    [~,Rtip,Ptip] = HTM_nume_mex(p,1,Lbody,params.r(n), ...
        params.K(n),params.A(n));
    Pbase = Pbase+Rbase*Ptip;
    Rbase = Rbase*Rtip;
end
bases(:,4) = Pbase;
end
