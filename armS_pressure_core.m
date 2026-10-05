function [M,C,G,dM] = armS_pressure_core(q,dq,params)
% Three-section recursive dynamics in six independent pressure coordinates.
% q = [p12;p13;p22;p23;p32;p33] (Pa); dq has the same order (Pa/s).
% Each section has uniform distributed mass mi(n)*dxi, xi in [0,1].
% Translational kinetic energy is included; rotational kinetic energy is not.
% G is the potential gradient for gravity acceleration params.g:
% V = -sum_n mi(n)*int(g.'*x_n(xi),xi=0..1), so M*qdd+C*dq+G=Q.
% params.L(:,n) = [base sensor offset; section length; tip sensor offset].
% Only params.L(2,n) enters mass integrals and section-to-section motion.
% The polynomial local rotation is treated as orthogonal in the local M block.

nq = 6;
assert(isequal(size(q),[nq,1]) && isequal(size(dq),[nq,1]));
assert(isequal(size(params.L),[3,3]));
assert(isequal(size(params.r),[3,1]) && isequal(size(params.K),[3,1]) ...
    && isequal(size(params.A),[3,1]) && isequal(size(params.mi),[3,1]) ...
    && isequal(size(params.g),[3,1]));
assert(all(params.L(2,:)>0) && all(params.r>0) ...
    && all(params.K>0) && all(params.A>0) && all(params.mi>0));

M = zeros(nq); dM = zeros(nq,nq,nq); G = zeros(nq,1);
R = eye(3);
Pq = zeros(3,nq); Rq = zeros(3,3,nq);
Pqq = zeros(3,nq,nq); Rqq = zeros(3,3,nq,nq);

for n = 1:3
    cur = 2*n-1:2*n;
    p = [0,q(cur).'];
    Lbody = [0;params.L(2,n);0];
    rn = params.r(n); Kn = params.K(n); An = params.A(n);
    m = params.mi(n);
    [mu,muq,muqq] = integratedPosition_nume(p,Lbody,rn,Kn,An);
    [S,Sq] = integratedPositionProduct_nume(p,Lbody,rn,Kn,An);
    [F,Fq] = integratedPositionDerivativeProduct_nume(p,Lbody,rn,Kn,An);
    [E,Eq] = integratedJacobianProduct_nume(p,Lbody,rn,Kn,An);

    % For an upstream coordinate j, J_j=Pq_j+Rq_j*p(xi).
    % For a local coordinate a, J_a=R*p_,a(xi).
    for i = 1:2*(n-1)
        Ai = Pq(:,i); Bi = Rq(:,:,i);
        G(i) = G(i)-m*(Ai+Bi*mu).'*params.g;
        for j = 1:2*(n-1)
            Aj = Pq(:,j); Bj = Rq(:,:,j);
            M(i,j) = M(i,j)+m*(Ai.'*Aj+Ai.'*Bj*mu ...
                +Aj.'*Bi*mu+trace(Bi.'*Bj*S));
            for h = 1:2*n
                Aih = Pqq(:,i,h); Bih = Rqq(:,:,i,h);
                Ajh = Pqq(:,j,h); Bjh = Rqq(:,:,j,h);
                muh = zeros(3,1); Sh = zeros(3);
                if h>2*(n-1)
                    a = h-2*(n-1);
                    muh = muq(:,a); Sh = Sq(:,:,a);
                end
                dM(i,j,h) = dM(i,j,h)+m*(Aih.'*Aj+Ai.'*Ajh ...
                    +Aih.'*Bj*mu+Ai.'*Bjh*mu+Ai.'*Bj*muh ...
                    +Ajh.'*Bi*mu+Aj.'*Bih*mu+Aj.'*Bi*muh ...
                    +trace((Bih.'*Bj+Bi.'*Bjh)*S+Bi.'*Bj*Sh));
            end
        end
        for a = 1:2
            j = cur(a);
            v = m*(Ai.'*R*muq(:,a)+trace(Bi.'*R*F(:,:,a)));
            M(i,j) = M(i,j)+v; M(j,i) = M(j,i)+v;
            for h = 1:2*n
                uq_h = zeros(3,1); Fh = zeros(3);
                if h>2*(n-1)
                    b = h-2*(n-1);
                    uq_h = muqq(:,a,b); Fh = Fq(:,:,a,b);
                end
                v = m*(Pqq(:,i,h).'*R*muq(:,a) ...
                    +Ai.'*Rq(:,:,h)*muq(:,a)+Ai.'*R*uq_h ...
                    +trace((Rqq(:,:,i,h).'*R+Bi.'*Rq(:,:,h))*F(:,:,a) ...
                    +Bi.'*R*Fh));
                dM(i,j,h) = dM(i,j,h)+v;
                dM(j,i,h) = dM(j,i,h)+v;
            end
        end
    end
    M(cur,cur) = M(cur,cur)+m*E;
    for a = 1:2
        dM(cur,cur,cur(a)) = dM(cur,cur,cur(a))+m*Eq(:,:,a);
        G(cur(a)) = G(cur(a))-m*(R*muq(:,a)).'*params.g;
    end

    % Propagate section-tip derivatives to the next section's base.
    [~,Rt,ptip] = HTM_nume_mex(p,1,Lbody,rn,Kn,An);
    [pj,rj,pjj,rjj] = LocalJacob_nume_mex(p,1,Lbody,rn,Kn,An);
    dp = zeros(3,nq); dr = zeros(3,3,nq);
    ddp = zeros(3,nq,nq); ddr = zeros(3,3,nq,nq);
    for a = 1:2
        rows = 3*(a-1)+(1:3);
        dp(:,cur(a)) = pj(:,a);
        dr(:,:,cur(a)) = rj(:,rows);
        for b = 1:2
            cols = 3*(b-1)+(1:3);
            ddp(:,cur(a),cur(b)) = pjj(rows,b);
            ddr(:,:,cur(a),cur(b)) = rjj(rows,cols);
        end
    end
    newPq = zeros(3,nq); newRq = zeros(3,3,nq);
    newPqq = zeros(3,nq,nq); newRqq = zeros(3,3,nq,nq);
    for i = 1:2*n
        newPq(:,i) = Pq(:,i)+Rq(:,:,i)*ptip+R*dp(:,i);
        newRq(:,:,i) = Rq(:,:,i)*Rt+R*dr(:,:,i);
        for h = 1:2*n
            newPqq(:,i,h) = Pqq(:,i,h)+Rqq(:,:,i,h)*ptip ...
                +Rq(:,:,i)*dp(:,h)+Rq(:,:,h)*dp(:,i)+R*ddp(:,i,h);
            newRqq(:,:,i,h) = Rqq(:,:,i,h)*Rt ...
                +Rq(:,:,i)*dr(:,:,h)+Rq(:,:,h)*dr(:,:,i)+R*ddr(:,:,i,h);
        end
    end
    R = R*Rt;
    Pq = newPq; Rq = newRq; Pqq = newPqq; Rqq = newRqq;
end

C = zeros(nq);
for i = 1:nq
    for j = 1:nq
        for h = 1:nq
            C(i,j) = C(i,j)+0.5*(dM(i,j,h)+dM(i,h,j) ...
                -dM(j,h,i))*dq(h);
        end
    end
end
end
