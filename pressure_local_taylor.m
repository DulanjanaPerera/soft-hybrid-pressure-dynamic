function [R,P,Pq,Rq,Pqq,Rqq] = pressure_local_taylor(p,xi,L,r,K,A)
%#codegen
% Compact pressure derivatives of the exported transform. Its rotation is
% exponential; its position uses the polynomial through W^12.
c = A*r/K;
W2 = [0,0,3*c/2;0,0,-sqrt(3)*c/2;-3*c/2,sqrt(3)*c/2,0];
W3 = [0,0,3*c/2;0,0,sqrt(3)*c/2;-3*c/2,-sqrt(3)*c/2,0];
W = p(2)*W2+p(3)*W3;
Wa = zeros(3,3,2);
Wa(:,:,1) = W2;
Wa(:,:,2) = W3;
B = xi*W;
Ba = xi*Wa;
R=zeros(3); P=zeros(3,1);
Pq=zeros(3,2); Rq=zeros(3,3,2);
Pqq=zeros(3,2,2); Rqq=zeros(3,3,2,2);
T=eye(3); Tq=zeros(3,3,2); Tqq=zeros(3,3,2,2);
ez=[0;0;1];
for order=0:30
    R=R+T;
    if order<=12
        P=P+L(2)*xi/(order+1)*(T*ez);
    end
    for a=1:2
        Rq(:,:,a)=Rq(:,:,a)+Tq(:,:,a);
        if order<=12
            Pq(:,a)=Pq(:,a)+L(2)*xi/(order+1)*(Tq(:,:,a)*ez);
        end
        for b=1:2
            Rqq(:,:,a,b)=Rqq(:,:,a,b)+Tqq(:,:,a,b);
            if order<=12
                Pqq(:,a,b)=Pqq(:,a,b)+L(2)*xi/(order+1)*(Tqq(:,:,a,b)*ez);
            end
        end
    end
    if order<30
        nextT=T*B/(order+1);
        nextTq=zeros(3,3,2); nextTqq=zeros(3,3,2,2);
        for a=1:2
            nextTq(:,:,a)=(Tq(:,:,a)*B+T*Ba(:,:,a))/(order+1);
            for b=1:2
                nextTqq(:,:,a,b)=(Tqq(:,:,a,b)*B ...
                    +Tq(:,:,a)*Ba(:,:,b)+Tq(:,:,b)*Ba(:,:,a))/(order+1);
            end
        end
        T=nextT; Tq=nextTq; Tqq=nextTqq;
    end
end
% Sensor extensions are excluded by the dynamics, which passes L(1)=L(3)=0.
end
