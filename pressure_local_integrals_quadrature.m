function [mu,muq,muqq,S,Sq,F,Fq,E,Eq] = ...
    pressure_local_integrals_quadrature(p,L,r,K,A)
%#codegen
% Integrate the same polynomial position/Jacobian/Hessian used by the
% symbolic exports. 16-point Gauss-Legendre is exact for polynomials through degree 31.
x = [0.0052995325041750307;0.0277124884633837;0.067184398806084122;0.1222977958224985;0.19106187779867811;0.27099161117138632;0.35919822461037054;0.45249374508118129;0.54750625491881877;0.64080177538962946;0.72900838882861363;0.80893812220132189;0.87770220417750155;0.93281560119391593;0.9722875115366163;0.99470046749582497];
w = [0.013576229705877019;0.031126761969323853;0.047579255841246296;0.062314485627767015;0.074797994408288382;0.08457825969750131;0.091301707522461806;0.094725305227534293;0.094725305227534293;0.091301707522461806;0.08457825969750131;0.074797994408288382;0.062314485627767015;0.047579255841246296;0.031126761969323853;0.013576229705877019];
mu=zeros(3,1); muq=zeros(3,2); muqq=zeros(3,2,2);
S=zeros(3); Sq=zeros(3,3,2);
F=zeros(3,3,2); Fq=zeros(3,3,2,2);
E=zeros(2); Eq=zeros(2,2,2);
for z=1:16
    [~,v,J,~,H,~]=pressure_local_taylor(p,x(z),L,r,K,A);
    wz=w(z);
    mu=mu+wz*v;
    muq=muq+wz*J;
    S=S+wz*(v*v.');
    E=E+wz*(J.'*J);
    for a=1:2
        ha=reshape(H(:,a,:),3,2);
        Sq(:,:,a)=Sq(:,:,a)+wz*(J(:,a)*v.'+v*J(:,a).');
        F(:,:,a)=F(:,:,a)+wz*(J(:,a)*v.');
        for b=1:2
            muqq(:,a,b)=muqq(:,a,b)+wz*ha(:,b);
            Fq(:,:,a,b)=Fq(:,:,a,b)+wz*(ha(:,b)*v.'+J(:,a)*J(:,b).');
            for c=1:2
                hc=reshape(H(:,c,:),3,2);
                Eq(a,b,c)=Eq(a,b,c)+wz*(hc(:,a).'*J(:,b)+J(:,a).'*hc(:,b));
            end
        end
    end
end
end
