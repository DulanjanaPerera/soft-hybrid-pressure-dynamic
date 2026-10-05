function fig = plot_armS_pressure_history(t,X)
% Plot pressure-equivalent shape coordinates in reusable Figure 2.
% These are not predicted measured chamber pressures in a passive run.
% t is N-by-1 (s), X is N-by-12; only X(:,1:6) is plotted.
assert(isvector(t) && ~isempty(t) && size(X,1)==numel(t) ...
    && size(X,2)>=6 && all(isfinite(X(:,1:6)),'all'));
fig = figure(2);
clf(fig);
set(fig,'Color','w','Position',[1040 100 850 700]);
layout = tiledlayout(fig,3,1,'TileSpacing','compact', ...
    'Padding','compact');
for n = 1:3
    ax = nexttile(layout,n);
    rows = 2*n-1:2*n;
    plot(ax,t(:),X(:,rows)/1000,'LineWidth',1.5);
    grid(ax,'on');
    ylabel(ax,sprintf('Module %d q (kPa)',n));
    legend(ax,sprintf('p_{%d2}',n),sprintf('p_{%d3}',n), ...
        'Location','best');
end
xlabel(layout,'Time (s)');
title(layout,'Pressure-equivalent shape coordinates by module');
end
