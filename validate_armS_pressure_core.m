function validate_armS_pressure_core()
% Validate three-section pressure dynamics against direct arm geometry.
% The direct reference differentiates products of local transforms by
% pressure perturbation and integrates their global point Jacobians.

params.L = [0.010,0.008,0.012; ...
            0.278,0.260,0.290; ...
            0.015,0.012,0.006];
params.r = [0.013;0.0125;0.014];
params.K = [0.0226;0.021;0.024];
params.A = pi*(params.r/2).^2;
params.mi = [0.10;0.08;0.12];
params.g = [0;0;-9.81];
dq = [100;-70;50;20;-40;80];
poses = [zeros(6,1), ...
    [4500;1500;3000;1000;2000;4000], ...
    [10000;0;2000;6000;5000;1000]];

for k = 1:size(poses,2)
    q = poses(:,k);
    [M,C,G,dM] = armS_pressure_core(q,dq,params);
    [Mq,Gq] = quadratureReference(q,params);
    dMfd = zeros(6,6,6);
    for h = 1:6
        step = zeros(6,1); step(h) = 20; % Pa
        Mp = armS_pressure_core(q+step,dq,params);
        Mm = armS_pressure_core(q-step,dq,params);
        dMfd(:,:,h) = (Mp-Mm)/40;
    end
    Mdot = zeros(6);
    for h = 1:6
        Mdot = Mdot+dM(:,:,h)*dq(h);
    end
    skewMatrix = Mdot-2*C;
    symmetry = norm(M-M.','fro')/max(norm(M,'fro'),eps);
    derivative = norm(dM(:)-dMfd(:))/max(norm(dMfd(:)),1e-20);
    skew = norm(skewMatrix+skewMatrix.','fro') ...
        /max(norm(Mdot,'fro')+2*norm(C,'fro'),1e-20);
    massIntegral = norm(M-Mq,'fro')/max(norm(Mq,'fro'),eps);
    gravityAbs = norm(G-Gq);
    rotationDefect = maxRotationDefect(q,params);
    [~,pd] = chol((M+M.')/2);
    fprintf(['Pose %d: symmetry %.3e, dM %.3e, skew %.3e, ' ...
        'quadrature M %.3e, gravity abs %.3e, rotation %.3e, chol %d\n'], ...
        k,symmetry,derivative,skew,massIntegral,gravityAbs, ...
        rotationDefect,pd);
    assert(all(isfinite(M(:))) && all(isfinite(C(:))) ...
        && all(isfinite(G(:))) && all(isfinite(dM(:))));
    assert(pd==0,'The distributed-mass matrix is not positive definite.');
    assert(symmetry<1e-10 && derivative<1e-4 && skew<1e-10, ...
        'Recursive mass derivative or Coriolis consistency failed.');
    assert(massIntegral<1e-5 && gravityAbs<1e-8, ...
        'Recursive dynamics disagree with direct global-point quadrature.');
    assert(rotationDefect<1e-5, ...
        'Polynomial rotation is too far from orthogonal for the E block.');

    if k == 2
        Gpotential = zeros(6,1);
        for h = 1:6
            step = zeros(6,1); step(h) = 20; % Pa
            Gpotential(h) = (potentialReference(q+step,params) ...
                -potentialReference(q-step,params))/40;
        end
        potentialError = norm(G-Gpotential) ...
            /max(norm(Gpotential),1e-12);
        fprintf('Pose %d: gravity vs potential gradient %.3e\n', ...
            k,potentialError);
        assert(potentialError<1e-5, ...
            'Gravity sign or potential gradient is inconsistent.');
    end

    % Changing a massless sensor offset must not change mechanical terms.
    offsetFree = params;
    offsetFree.L([1,3],:) = 0;
    [M0,C0,G0,dM0] = armS_pressure_core(q,dq,offsetFree);
    assert(isequal(M,M0) && isequal(C,C0) ...
        && isequal(G,G0) && isequal(dM,dM0));
end
fprintf('All three-section pressure-core checks passed.\n');
end

function defect = maxRotationDefect(q,params)
T = eye(4);
defect = 0;
for n = 1:3
    cur = 2*n-1:2*n;
    Lbody = [0;params.L(2,n);0];
    T = T*HTM_nume_mex([0,q(cur).'],1,Lbody, ...
        params.r(n),params.K(n),params.A(n));
    R = T(1:3,1:3);
    defect = max(defect,norm(R.'*R-eye(3),'fro'));
end
end

function V = potentialReference(q,params)
% Direct gravitational potential of the uniformly distributed backbone.
n = 20;
k = (1:n-1).';
b = k./sqrt(4*k.^2-1);
[U,D] = eig(diag(b,1)+diag(b,-1));
xi = (diag(D)+1)/2;
w = (U(1,:).^2).';
V = 0;
for section = 1:3
    for s = 1:n
        x = globalPoint(q,section,xi(s),params);
        V = V-params.mi(section)*w(s)*(params.g.'*x);
    end
end
end

function [M,G] = quadratureReference(q,params)
% Independent global point Jacobians from full transform products.
n = 20;
k = (1:n-1).';
b = k./sqrt(4*k.^2-1);
[V,D] = eig(diag(b,1)+diag(b,-1));
xi = (diag(D)+1)/2;
w = (V(1,:).^2).';
M = zeros(6); G = zeros(6,1);
for section = 1:3
    for s = 1:n
        J = zeros(3,6);
        for h = 1:6
            step = zeros(6,1); step(h) = 2; % Pa
            xp = globalPoint(q+step,section,xi(s),params);
            xm = globalPoint(q-step,section,xi(s),params);
            J(:,h) = (xp-xm)/4;
        end
        M = M+params.mi(section)*w(s)*(J.'*J);
        G = G-params.mi(section)*w(s)*(J.'*params.g);
    end
end
end

function x = globalPoint(q,section,xi,params)
% Direct composition of section transforms, without recursive derivatives.
T = eye(4);
for n = 1:section
    cur = 2*n-1:2*n;
    p = [0,q(cur).'];
    Lbody = [0;params.L(2,n);0];
    if n == section
        localXi = xi;
    else
        localXi = 1;
    end
    T = T*HTM_nume_mex(p,localXi,Lbody, ...
        params.r(n),params.K(n),params.A(n));
end
x = T(1:3,4);
end
