function validate_armS_pressure_dynamics()
% Check elastic-law provenance, potential gradient, and passive release.
[t,X,params,d] = simulate_armS_pressure_passive();
[~,~,H] = pressure_elastic_force(zeros(6,1),params);
assert(norm(H-H.','fro')<1e-20 && min(eig(H))>0);

poses = [X(1,1:6).',X(round(end/2),1:6).',X(end,1:6).'];
maxElasticError = 0; maxPotentialError = 0;
for j = 1:size(poses,2)
    q = poses(:,j);
    Fel = pressure_elastic_force(q,params);
    % The third output from the core is the gravity potential gradient.
    [~,~,G] = armS_pressure_core(q,zeros(6,1),params);
    potentialFD = zeros(6,1);
    for h = 1:6
        step = zeros(6,1); step(h) = 2; % Pa
        [~,~,Up,Vp] = armS_pressure_energy(q+step,zeros(6,1),params);
        [~,~,Um,Vm] = armS_pressure_energy(q-step,zeros(6,1),params);
        potentialFD(h) = ((Up+Vp)-(Um+Vm))/4;
    end
    potentialError = norm(potentialFD-(Fel+G)) ...
        /max(norm(Fel+G),1e-12);
    maxPotentialError = max(maxPotentialError,potentialError);
    for n = 1:3
        rows = 2*n-1:2*n;
        p = [0;q(rows)];
        original = G_and_K_pressure_2DoF_MSF(p,params.mi(n), ...
            [0;params.L(2,n);0],params.r(n),params.K(n), ...
            params.A(n),params.k(:,n),0);
        elasticError = norm(original-Fel(rows)) ...
            /max(norm(Fel(rows)),1e-12);
        maxElasticError = max(maxElasticError,elasticError);
    end
end
energyDrop = d.energy(1)-d.energy(end);
balanceError = abs(d.energyBalanceResidual)/max(energyDrop,1e-12);
fprintf(['samples %d, max |q| %.1f Pa, elastic match %.3e, ' ...
    'potential gradient %.3e, E0 %.9g J, Ef %.9g J, ' ...
    'balance %.3e, max sampled rise %.3e J\n'], ...
    numel(t),max(abs(X(:,1:6)),[],'all'),maxElasticError, ...
    maxPotentialError,d.energy(1),d.energy(end),balanceError, ...
    max(diff(d.energy)));
assert(all(isfinite(X(:))) && all(isfinite(d.energy)));
assert(maxElasticError<1e-7 && maxPotentialError<1e-4);
assert(all(diff(d.energy)<=1e-8) && energyDrop>0);
assert(balanceError<5e-3);
end
