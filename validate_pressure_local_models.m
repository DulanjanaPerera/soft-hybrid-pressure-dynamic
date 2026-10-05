function validate_pressure_local_models()
% Independent local pressure-derivative and section-integral checks.
% Pressure order is [0,p2,p3] in Pa. The integration measure is dxi.

validate_pressure_local_interfaces();
r = 0.013;
K = 0.0226;
A = pi*(r/2)^2;
pressures = {[0,0,0], [0,4500,1500], [0,2000,8000], ...
    [0,10000,0]};
lengths = {[0.01;0.278;0.015], [0.008;0.260;0.012], ...
    [0.012;0.290;0.006], [0.01;0.278;0.015]};
xis = [0.35, 0.67, 1.0, 0.5];

for k = 1:numel(pressures)
    p = pressures{k};
    L = lengths{k};
    xi = xis(k);
    [~,R,P] = HTM_nume_mex(p,xi,L,r,K,A);
    [J,RJ,H,RH] = LocalJacob_nume_mex(p,xi,L,r,K,A);
    if k == 1
        assert(norm(R-eye(3),'fro') < 1e-12);
        assert(norm(P-[0;0;L(1)+xi*L(2)+L(3)]) < 1e-12);
    end
    firstP = 0;
    firstR = 0;
    secondP = 0;
    secondR = 0;
    for b = 1:2
        h = 2; % Pa; pressure-scale step for centered differences.
        pp = p; pm = p;
        pp(b+1) = pp(b+1)+h;
        pm(b+1) = pm(b+1)-h;
        [~,Rp,Pp] = HTM_nume_mex(pp,xi,L,r,K,A);
        [~,Rm,Pm] = HTM_nume_mex(pm,xi,L,r,K,A);
        firstP = max(firstP,scaledError(J(:,b),(Pp-Pm)/(2*h),1e-12));
        cols = 3*(b-1)+(1:3);
        firstR = max(firstR,scaledError(RJ(:,cols),(Rp-Rm)/(2*h),1e-12));
        h2 = 20; % Pa; larger step reduces cancellation in second derivatives.
        pp2 = p; pm2 = p;
        pp2(b+1) = pp2(b+1)+h2;
        pm2(b+1) = pm2(b+1)-h2;
        [Jp,RJp] = LocalJacob_nume_mex(pp2,xi,L,r,K,A);
        [Jm,RJm] = LocalJacob_nume_mex(pm2,xi,L,r,K,A);
        for a = 1:2
            rows = 3*(a-1)+(1:3);
            secondP = max(secondP,scaledError(H(rows,b), ...
                (Jp(:,a)-Jm(:,a))/(2*h2),1e-12));
            cols = 3*(b-1)+(1:3);
            aCols = 3*(a-1)+(1:3);
            secondR = max(secondR,scaledError(RH(rows,cols), ...
                (RJp(:,aCols)-RJm(:,aCols))/(2*h2),1e-12));
        end
    end
    fprintf('Pose %d: pressure FD P %.3e, R %.3e, Ppp %.3e, Rpp %.3e\n', ...
        k,firstP,firstR,secondP,secondR);

    % Sensor offsets are geometry only. Distributed mass uses the backbone.
    Lbody = [0;L(2);0];
    [mu,muq,muqq] = integratedPosition_nume(p,Lbody,r,K,A);
    [S,Sq] = integratedPositionProduct_nume(p,Lbody,r,K,A);
    [F,Fq] = integratedPositionDerivativeProduct_nume(p,Lbody,r,K,A);
    [E,Eq] = integratedJacobianProduct_nume(p,Lbody,r,K,A);
    [muRef,muqRef,muqqRef,SRef,SqRef,FRef,FqRef,ERef,EqRef] = ...
        quadratureReference(p,Lbody,r,K,A);
    analytic = {mu,muq,muqq,S,Sq,F,Fq,E,Eq};
    reference = {muRef,muqRef,muqqRef,SRef,SqRef,FRef,FqRef,ERef,EqRef};
    names = {'mu','muq','muqq','S','Sq','F','Fq','E','Eq'};
    errors = zeros(1,numel(names));
    for i = 1:numel(names)
        errors(i) = scaledError(analytic{i},reference{i},1e-12);
    end
    identity = 0;
    for a = 1:2
        identity = max(identity,scaledError(Sq(:,:,a), ...
            F(:,:,a)+F(:,:,a).',1e-12));
    end
    fprintf('Pose %d: quadrature errors ',k);
    for i = 1:numel(names)
        fprintf('%s %.3e ',names{i},errors(i));
    end
    fprintf('| Sq=F+F'' %.3e\n',identity);
    assert(max([firstP firstR secondP secondR]) < 1e-5, ...
        'Pressure derivatives disagree with centered differences.');
    assert(max([errors identity]) < 1e-6, ...
        'Analytic moments disagree with direct quadrature.');
    assert(all(isfinite(R(:))) && all(isfinite(P(:))));
end
fprintf('All local pressure-model checks passed.\n');
end

function [mu,muq,muqq,S,Sq,F,Fq,E,Eq] = quadratureReference(p,L,r,K,A)
% 24-point Gauss-Legendre integration of the point geometry and derivatives.
n = 24;
k = (1:n-1).';
b = k./sqrt(4*k.^2-1);
[V,D] = eig(diag(b,1)+diag(b,-1));
xi = (diag(D)+1)/2;
w = (V(1,:).^2).';
mu = zeros(3,1); muq = zeros(3,2); muqq = zeros(3,2,2);
S = zeros(3); Sq = zeros(3,3,2);
F = zeros(3,3,2); Fq = zeros(3,3,2,2);
E = zeros(2); Eq = zeros(2,2,2);
for s = 1:n
    [~,~,P] = HTM_nume_mex(p,xi(s),L,r,K,A);
    [J,~,H] = LocalJacob_nume_mex(p,xi(s),L,r,K,A);
    H3 = zeros(3,2,2);
    for a = 1:2
        rows = 3*(a-1)+(1:3);
        for c = 1:2
            H3(:,a,c) = H(rows,c);
        end
    end
    mu = mu+w(s)*P;
    muq = muq+w(s)*J;
    muqq = muqq+w(s)*H3;
    S = S+w(s)*(P*P.');
    E = E+w(s)*(J.'*J);
    for a = 1:2
        F(:,:,a) = F(:,:,a)+w(s)*(J(:,a)*P.');
        Sq(:,:,a) = Sq(:,:,a)+w(s)*(J(:,a)*P.'+P*J(:,a).');
        for c = 1:2
            Fq(:,:,a,c) = Fq(:,:,a,c)+w(s)*(H3(:,a,c)*P.' ...
                +J(:,a)*J(:,c).');
        end
        Eq(:,:,a) = Eq(:,:,a)+w(s)*(H3(:,:,a).'*J+J.'*H3(:,:,a));
    end
end
end

function e = scaledError(actual,expected,scaleFloor)
e = norm(actual(:)-expected(:))/max(norm(expected(:)),scaleFloor);
end
