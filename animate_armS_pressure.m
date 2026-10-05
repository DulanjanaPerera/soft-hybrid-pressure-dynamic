function outputFile = animate_armS_pressure(t,X,params,opts)
% Animate three connected pressure-coordinate sections as colored tubes.
% Optional opts: saveVideo, fileName, frameStride, fps, nXi, nCircle,
% bodyRadius, visible, closeFigure. X(:,1:6) contains the pressures (Pa).
if nargin<4, opts = struct(); end
assert(isvector(t) && size(X,1)==numel(t) && size(X,2)>=6);
assert(all(isfinite(t(:))) && all(isfinite(X(:,1:6)),'all'));
saveVideo = option(opts,'saveVideo',false);
frameStride = option(opts,'frameStride',8);
fps = option(opts,'fps',25);
nXi = option(opts,'nXi',31);
nCircle = option(opts,'nCircle',16);
radius = option(opts,'bodyRadius',0.015);
visible = option(opts,'visible','on');
closeFigure = option(opts,'closeFigure',false);
outputFile = option(opts,'fileName',fullfile('results', ...
    'armS_pressure_passive.mp4'));
assert(frameStride>=1 && frameStride==round(frameStride) && ...
    nXi>=2 && nXi==round(nXi) && nCircle>=8 && ...
    nCircle==round(nCircle) && fps>0 && radius>0);
if saveVideo && strcmpi(visible,'off')
    error('Video capture requires a visible MATLAB figure.');
end

frames = 1:frameStride:numel(t);
if frames(end)~=numel(t), frames(end+1) = numel(t); end
xi = linspace(0,1,nXi);
colors = [0.18 0.48 0.85;0.13 0.65 0.47;0.92 0.51 0.18];
fig = figure('Color','w','Visible',visible,'Position',[100 100 900 700]);
ax = axes('Parent',fig); hold(ax,'on'); grid(ax,'on');
axis(ax,'equal'); view(ax,35,22);
xlabel(ax,'X (m)'); ylabel(ax,'Y (m)'); zlabel(ax,'Z (m)');
mins = zeros(3,1); maxs = zeros(3,1);
for k = frames
    Psample = armS_pressure_geometry(X(k,1:6).',params,xi);
    points = reshape(permute(Psample,[2,1,3]),3,[]);
    mins = min(mins,min(points,[],2));
    maxs = max(maxs,max(points,[],2));
end
center = (mins+maxs)/2;
halfSpan = max([(maxs-mins).'/2+3*radius,0.15]);
xlim(ax,center(1)+[-halfSpan halfSpan]);
ylim(ax,center(2)+[-halfSpan halfSpan]);
zlim(ax,center(3)+[-halfSpan halfSpan]);
surfaces = gobjects(3,1);
for n = 1:3
    surfaces(n) = surf(ax,nan(nXi,nCircle+1),nan(nXi,nCircle+1), ...
        nan(nXi,nCircle+1),'FaceColor',colors(n,:), ...
        'EdgeColor','none','FaceLighting','gouraud');
end
tip = plot3(ax,0,0,0,'ko','MarkerFaceColor','k','MarkerSize',6);
camlight(ax,'headlight');
if saveVideo
    [folder,~,ext] = fileparts(outputFile);
    if isempty(ext), outputFile = [outputFile,'.mp4']; end
    if ~isempty(folder) && ~isfolder(folder), mkdir(folder); end
    try
        video = VideoWriter(outputFile,'MPEG-4');
    catch
        warning('MPEG-4 unavailable; using Motion JPEG AVI.');
        [folder,name] = fileparts(outputFile);
        outputFile = fullfile(folder,[name,'.avi']);
        video = VideoWriter(outputFile,'Motion JPEG AVI');
    end
    video.FrameRate = fps;
    open(video);
end
try
    for k = frames
        [P,R,bases] = armS_pressure_geometry(X(k,1:6).',params,xi);
        for n = 1:3
            [xx,yy,zz] = tube(P(:,:,n),R(:,:,:,n),radius,nCircle);
            set(surfaces(n),'XData',xx,'YData',yy,'ZData',zz);
        end
        set(tip,'XData',bases(1,4),'YData',bases(2,4), ...
            'ZData',bases(3,4));
        title(ax,sprintf('Three-section pressure dynamics | t = %.3f s',t(k)));
        drawnow;
        if saveVideo, writeVideo(video,getframe(fig)); end
    end
catch exception
    if saveVideo, close(video); end
    rethrow(exception);
end
if saveVideo, close(video); end
if closeFigure, close(fig); end
end

function [xx,yy,zz] = tube(P,R,radius,nCircle)
theta = linspace(0,2*pi,nCircle+1);
xx = zeros(size(P,1),nCircle+1);
yy = xx; zz = xx;
for j = 1:size(P,1)
    ring = P(j,:).' + radius*(R(:,1,j)*cos(theta) ...
        +R(:,2,j)*sin(theta));
    xx(j,:) = ring(1,:);
    yy(j,:) = ring(2,:);
    zz(j,:) = ring(3,:);
end
end

function value = option(opts,name,defaultValue)
if isfield(opts,name), value = opts.(name);
else, value = defaultValue;
end
end
