function validate_armS_pressure_geometry()
% Check connected physical backbones and sensor-offset independence.
[~,X,params] = simulate_armS_pressure_passive();
xi = linspace(0,1,11);
poses = [zeros(6,1),X(1,1:6).',X(round(end/2),1:6).', ...
    X(end,1:6).'];
maxJointGap = 0; maxSensorEffect = 0;
for k = 1:size(poses,2)
    q = poses(:,k);
    [P,R,bases] = armS_pressure_geometry(q,params,xi);
    for n = 1:3
        gap = norm(P(end,:,n).'-bases(:,n+1));
        maxJointGap = max(maxJointGap,gap);
        if n<3
            gap = norm(P(end,:,n)-P(1,:,n+1));
            maxJointGap = max(maxJointGap,gap);
        end
    end
    altered = params;
    altered.L([1,3],:) = altered.L([1,3],:)+0.1;
    [Pa,Ra,Ba] = armS_pressure_geometry(q,altered,xi);
    maxSensorEffect = max(maxSensorEffect,max(abs(Pa(:)-P(:))));
    maxSensorEffect = max(maxSensorEffect,max(abs(Ra(:)-R(:))));
    maxSensorEffect = max(maxSensorEffect,max(abs(Ba(:)-bases(:))));
    assert(all(isfinite(P(:))) && all(isfinite(R(:))));
    if k==1
        expected = cumsum([0,params.L(2,:)]);
        assert(norm(bases(3,:)-expected)<1e-12);
        assert(max(abs(bases(1:2,:)),[],'all')<1e-12);
    end
end
fprintf('geometry: max joint gap %.3e m, sensor-offset effect %.3e\n', ...
    maxJointGap,maxSensorEffect);
assert(maxJointGap<1e-10 && maxSensorEffect==0);
end
