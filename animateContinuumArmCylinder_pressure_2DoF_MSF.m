function outputFile = animateContinuumArmCylinder_pressure_2DoF_MSF(t, X, params, opts)
% animateContinuumArmCylinder_pressure_2DoF_MSF
%
% Draws the soft continuum arm as a 3D cylindrical/tubular body and,
% optionally, records the animation to a video file.
%
% Required inputs:
%   t      : time vector from ode15s                  [N x 1]
%   X      : state history from ode15s                [N x 4]
%            X(:,1:2) are the two independent pressure-space coordinates
%   params : same parameter struct used in simulation
%
% Optional input:
%   opts   : visualization/recording options
%
% This function assumes the following function is available on the path:
%   backbonePos_pressure_xi_2DoF_MSF(p, params.L, params.r, xi, params.K, params.A)

    if nargin < 4
        opts = struct();
    end

    % -------------------- Options --------------------
    outputFile     = getOpt(opts, 'fileName',      'soft_arm_cylinder_animation.mp4');
    saveVideo      = getOpt(opts, 'saveVideo',     true);

    bodyRadius     = getOpt(opts, 'bodyRadius',    0.015);
    nXi            = getOpt(opts, 'nXi',           60);
    nCircle        = getOpt(opts, 'nCircle',       24);
    frameStride    = getOpt(opts, 'frameStride',   1);
    fps            = getOpt(opts, 'fps',           30);
    pauseTime      = getOpt(opts, 'pauseTime',     0.01);
    showCenterline = getOpt(opts, 'showCenterline', false);

    bodyColor      = getOpt(opts, 'bodyColor',     [0.20 0.55 0.90]);
    edgeColor      = getOpt(opts, 'edgeColor',     [0.00 0.00 0.00]);
    edgeAlpha      = getOpt(opts, 'edgeAlpha',     0.65);
    faceAlpha      = getOpt(opts, 'faceAlpha',     1.00);
    edgeLineWidth  = getOpt(opts, 'edgeLineWidth', 0.35);

    xLimVal        = getOpt(opts, 'xlim',          [-0.3 0.3]);
    yLimVal        = getOpt(opts, 'ylim',          [-0.3 0.3]);
    zLimVal        = getOpt(opts, 'zlim',          [-0.1 0.3]);
    viewAngle      = getOpt(opts, 'view',          [45 25]);

    % -------------------- Basic checks --------------------
    nFrames = min(length(t), size(X,1));

    if size(X,2) < 2
        error('X must contain at least two state columns: X(:,1:2).');
    end

    frameStride = max(1, round(frameStride));
    nXi         = max(3, round(nXi));
    nCircle     = max(8, round(nCircle));

    xi = linspace(0, 1, nXi);

    % -------------------- Figure setup --------------------
    fig = figure('Color', 'w');
    clf(fig);

    ax = axes('Parent', fig);
    hold(ax, 'on');
    grid(ax, 'on');
    axis(ax, 'equal');

    xlabel(ax, 'X [m]');
    ylabel(ax, 'Y [m]');
    zlabel(ax, 'Z [m]');
    title(ax, 'Soft continuum arm cylindrical animation');

    xlim(ax, xLimVal);
    ylim(ax, yLimVal);
    zlim(ax, zLimVal);
    view(ax, viewAngle);
    rotate3d(ax, 'on');

    % -------------------- Initial geometry --------------------
    P0 = computeBackbonePoints(X(1,:), params, xi);

    [Xs, Ys, Zs, ringBase, ringTip] = makeTubeMesh(P0, bodyRadius, nCircle);
    [baseVertices, baseFaces] = makeCapPatch(P0(1,:), ringBase, true);
    [tipVertices,  tipFaces]  = makeCapPatch(P0(end,:), ringTip, false);

    % Main cylindrical surface
    hBody = surf(ax, Xs, Ys, Zs, ...
        'FaceColor', bodyColor, ...
        'EdgeColor', edgeColor, ...
        'EdgeAlpha', edgeAlpha, ...
        'LineWidth', edgeLineWidth, ...
        'FaceAlpha', faceAlpha, ...
        'FaceLighting', 'gouraud', ...
        'AmbientStrength', 0.35, ...
        'DiffuseStrength', 0.75, ...
        'SpecularStrength', 0.20);

    % Base cap
    hBaseCap = patch(ax, ...
        'Vertices', baseVertices, ...
        'Faces', baseFaces, ...
        'FaceColor', bodyColor, ...
        'EdgeColor', edgeColor, ...
        'EdgeAlpha', edgeAlpha, ...
        'LineWidth', edgeLineWidth, ...
        'FaceAlpha', faceAlpha, ...
        'FaceLighting', 'gouraud', ...
        'AmbientStrength', 0.35, ...
        'DiffuseStrength', 0.75, ...
        'SpecularStrength', 0.20);

    % Tip cap
    hTipCap = patch(ax, ...
        'Vertices', tipVertices, ...
        'Faces', tipFaces, ...
        'FaceColor', bodyColor, ...
        'EdgeColor', edgeColor, ...
        'EdgeAlpha', edgeAlpha, ...
        'LineWidth', edgeLineWidth, ...
        'FaceAlpha', faceAlpha, ...
        'FaceLighting', 'gouraud', ...
        'AmbientStrength', 0.35, ...
        'DiffuseStrength', 0.75, ...
        'SpecularStrength', 0.20);

    if showCenterline
        hCenter = plot3(ax, P0(:,1), P0(:,2), P0(:,3), ...
            'k-', 'LineWidth', 1.0);
    else
        hCenter = [];
    end

    hTip = plot3(ax, P0(end,1), P0(end,2), P0(end,3), ...
        'o', ...
        'MarkerSize', 7, ...
        'MarkerFaceColor', 'r', ...
        'MarkerEdgeColor', 'k');

    % -------------------- Lighting --------------------
    lighting(ax, 'gouraud');
    material(ax, 'dull');

    camlight(ax, 'headlight');
    camlight(ax, 'right');
    camlight(ax, 'left');

    light(ax, ...
        'Position', [0.5 -0.8 0.8], ...
        'Style', 'infinite');

    light(ax, ...
        'Position', [-0.5 0.8 0.6], ...
        'Style', 'infinite');

    % -------------------- Video setup --------------------
    if saveVideo
        [folderName, ~, ext] = fileparts(outputFile);

        if isempty(ext)
            outputFile = [outputFile, '.mp4'];
        end

        if ~isempty(folderName) && ~isfolder(folderName)
            mkdir(folderName);
        end

        try
            vid = VideoWriter(outputFile, 'MPEG-4');
        catch
            warning('MPEG-4 is not available. Saving as Motion JPEG AVI instead.');
            outputFile = replaceFileExtension(outputFile, '.avi');
            vid = VideoWriter(outputFile, 'Motion JPEG AVI');
        end

        vid.FrameRate = fps;
        open(vid);
    end

    % -------------------- Animation loop --------------------
    for k = 1:frameStride:nFrames

        P = computeBackbonePoints(X(k,:), params, xi);

        [Xs, Ys, Zs, ringBase, ringTip] = makeTubeMesh(P, bodyRadius, nCircle);
        [baseVertices, ~] = makeCapPatch(P(1,:), ringBase, true);
        [tipVertices,  ~] = makeCapPatch(P(end,:), ringTip, false);

        set(hBody, ...
            'XData', Xs, ...
            'YData', Ys, ...
            'ZData', Zs);

        set(hBaseCap, ...
            'Vertices', baseVertices);

        set(hTipCap, ...
            'Vertices', tipVertices);

        if showCenterline
            set(hCenter, ...
                'XData', P(:,1), ...
                'YData', P(:,2), ...
                'ZData', P(:,3));
        end

        set(hTip, ...
            'XData', P(end,1), ...
            'YData', P(end,2), ...
            'ZData', P(end,3));

        title(ax, sprintf('Soft continuum arm cylinder | frame %d / %d | t = %.3f s', ...
            k, nFrames, t(k)));

        drawnow;

        if saveVideo
            writeVideo(vid, getframe(fig));
        end

        if pauseTime > 0
            pause(pauseTime);
        end
    end

    if saveVideo
        close(vid);
        fprintf('Animation saved to: %s\n', outputFile);
    end
end

% ========================================================================
% Local helper functions
% ========================================================================

function P = computeBackbonePoints(Xrow, params, xi)
% Computes backbone centerline points for one simulation frame.

    p = zeros(3,1);
    p(2:3) = Xrow(1:2).';

    P = zeros(length(xi), 3);

    for j = 1:length(xi)
        pos = backbonePos_pressure_xi_2DoF_MSF( ...
            p, params.L, params.r, xi(j), params.K, params.A);

        P(j,:) = pos(:).';
    end
end

function [Xs, Ys, Zs, ringBase, ringTip] = makeTubeMesh(P, radius, nCircle)
% Creates a tube mesh around the 3D centerline P.
%
% P      : [nXi x 3] centerline
% radius : tube radius

    n = size(P,1);
    theta = linspace(0, 2*pi, nCircle + 1);

    T = computeTangents(P);
    [N, B] = computeNormalFrames(T);

    Xs = zeros(n, nCircle + 1);
    Ys = zeros(n, nCircle + 1);
    Zs = zeros(n, nCircle + 1);

    ringBase = zeros(nCircle, 3);
    ringTip  = zeros(nCircle, 3);

    for i = 1:n
        for c = 1:(nCircle + 1)
            radialVector = radius * (cos(theta(c))*N(i,:) + sin(theta(c))*B(i,:));
            point = P(i,:) + radialVector;

            Xs(i,c) = point(1);
            Ys(i,c) = point(2);
            Zs(i,c) = point(3);

            if i == 1 && c <= nCircle
                ringBase(c,:) = point;
            end

            if i == n && c <= nCircle
                ringTip(c,:) = point;
            end
        end
    end
end

function [vertices, faces] = makeCapPatch(centerPoint, ringPoints, flipNormal)
% Creates triangular fan cap geometry.
%
% vertices:
%   row 1      : cap center
%   rows 2:end : circular ring points
%
% faces:
%   triangular fan faces

    nCircle = size(ringPoints, 1);

    vertices = [centerPoint; ringPoints];

    faces = zeros(nCircle, 3);

    for c = 1:nCircle
        cNext = c + 1;
        if cNext > nCircle
            cNext = 1;
        end

        if flipNormal
            faces(c,:) = [1, cNext + 1, c + 1];
        else
            faces(c,:) = [1, c + 1, cNext + 1];
        end
    end
end

function T = computeTangents(P)
% Computes approximate tangent vectors along a 3D curve.

    n = size(P,1);
    T = zeros(n,3);

    for i = 1:n
        if i == 1
            dP = P(2,:) - P(1,:);
        elseif i == n
            dP = P(n,:) - P(n-1,:);
        else
            dP = P(i+1,:) - P(i-1,:);
        end

        T(i,:) = safeNormalize(dP, [0 0 1]);
    end
end

function [N, B] = computeNormalFrames(T)
% Computes smooth local normal/binormal frames along the centerline.

    n = size(T,1);
    N = zeros(n,3);
    B = zeros(n,3);

    ref = [0 0 1];

    if abs(dot(ref, T(1,:))) > 0.90
        ref = [0 1 0];
    end

    N(1,:) = safeNormalize(cross(ref, T(1,:)), [1 0 0]);
    B(1,:) = safeNormalize(cross(T(1,:), N(1,:)), [0 1 0]);

    for i = 2:n
        Ni = N(i-1,:) - dot(N(i-1,:), T(i,:)) * T(i,:);

        if norm(Ni) < 1e-10
            ref = [0 0 1];
            if abs(dot(ref, T(i,:))) > 0.90
                ref = [0 1 0];
            end
            Ni = cross(ref, T(i,:));
        end

        N(i,:) = safeNormalize(Ni, N(i-1,:));
        B(i,:) = safeNormalize(cross(T(i,:), N(i,:)), B(i-1,:));
    end
end

function v = safeNormalize(v, fallback)
% Normalizes a vector and uses a fallback if the vector is nearly zero.

    nv = norm(v);

    if nv < 1e-12 || any(~isfinite(v))
        v = fallback;
    else
        v = v ./ nv;
    end
end

function val = getOpt(opts, fieldName, defaultVal)
% Reads an option from opts if it exists; otherwise returns defaultVal.

    if isfield(opts, fieldName)
        val = opts.(fieldName);
    else
        val = defaultVal;
    end
end

function fileNameOut = replaceFileExtension(fileNameIn, newExt)
% Replaces the extension of a file name.

    [folderName, baseName, ~] = fileparts(fileNameIn);

    if isempty(folderName)
        fileNameOut = [baseName, newExt];
    else
        fileNameOut = fullfile(folderName, [baseName, newExt]);
    end
end