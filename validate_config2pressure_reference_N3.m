function validate_config2pressure_reference_N3()
% Round-trip and derivative checks for the section-angle reference map.
params.r = [0.013;0.014;0.012];
params.A = pi*(params.r/2).^2;
params.K = [0.01443936;0.017;0.013];
theta = [0,0.2,-0.4;0.7,-1.0,1.2;-0.3,0.5,-1.5];
phi = [0,0.4,0.8;0.2,0.5,0.9;0.1,0.3,0.7];
q = config2pressure_reference_N3(theta,phi,params);
worstAngle = 0;
for j=1:size(theta,2)
    for n=1:3
        rows=2*n-1:2*n;
        c=pressure2config([0;q(rows,j)], ...
            params.r(n),params.K(n),params.A(n));
        assert(abs(c(2)-phi(n,j))<1e-12);
        if phi(n,j)>0
            angleError=atan2(sin(c(1)-theta(n,j)), ...
                cos(c(1)-theta(n,j)));
            worstAngle=max(worstAngle,abs(angleError));
        end
    end
end
assert(worstAngle<1e-12);

theta0=[0.1;-0.5;0.8]; phi0=[0.3;0.4;0.5];
thetaDot=[0.4;-0.2;0.3]; phiDot=[0.1;0.2;-0.1];
thetaDDot=[-0.2;0.3;0.1]; phiDDot=[0.05;-0.1;0.2];
[q,dq,ddq]=config2pressure_reference_N3(theta0,phi0,params, ...
    thetaDot,phiDot,thetaDDot,phiDDot);
h=1e-4;
qp=config2pressure_reference_N3( ...
    theta0+h*thetaDot+0.5*h^2*thetaDDot, ...
    phi0+h*phiDot+0.5*h^2*phiDDot,params);
qm=config2pressure_reference_N3( ...
    theta0-h*thetaDot+0.5*h^2*thetaDDot, ...
    phi0-h*phiDot+0.5*h^2*phiDDot,params);
rateError=max(abs((qp-qm)/(2*h)-dq));
accelError=max(abs((qp-2*q+qm)/h^2-ddq));
assert(rateError<1e-4 && accelError<1e-3);
fprintf('Reference map passed: angle %.3e rad, rate %.3e Pa/s, acceleration %.3e Pa/s^2.\n', ...
    worstAngle,rateError,accelError);
end
