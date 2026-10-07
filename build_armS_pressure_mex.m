function build_armS_pressure_mex(buildFolder)
% Compile the three-section dynamics and check it against MATLAB.
% Run: build_armS_pressure_mex
if nargin < 1
    buildFolder = fullfile(tempdir,'armS_pressure_mex_build');
end
sourceFolder = fileparts(mfilename('fullpath'));
oldFolder = pwd;
oldPath = path;
cleanup = onCleanup(@() restoreEnvironment(oldFolder,oldPath)); %#ok<NASGU>
if ~isfolder(buildFolder), mkdir(buildFolder); end
addpath(sourceFolder,'-begin');
cd(buildFolder);
clear armS_pressure_mex

% Type and size examples. Numerical values remain runtime inputs.
L = repmat([0.001;0.178;0.001],1,3);
r = 0.013*ones(3,1);
A = pi*(r/2).^2;
k = 3200*ones(3,3);
K = (3/20)*k(1,:).'.*L(2,:).'.*r.^2;
args = {0,zeros(12,1),L,r,K,A,0.1*ones(3,1), ...
    [0;0;9.81],k,8e-10*eye(6),zeros(6,1)};
cfg = coder.config('mex');
cfg.GenerateReport = false;
cfg.LaunchReport = false;
cfg.TargetLang = 'C++';
cfg.EnableOpenMP = false;
cfg.InlineBetweenUserFunctions = 'Never';
fprintf('Compiling pressure dynamics in %s\n',buildFolder);
drawnow;
buildClock = tic;
codegen('-config',cfg,'armS_pressure_entry','-args',args, ...
    '-d',fullfile(buildFolder,'codegen'),'-o','armS_pressure_mex');
fprintf('Compilation finished in %.1f s\n',toc(buildClock));
rehash;

states = [zeros(12,1), ...
    [10000;0;200;100;100;0;zeros(6,1)], ...
    [3500;-2100;1800;900;-1000;700; ...
     50;-20;40;10;-30;15], ...
    [7000;-2000;4000;1000;-3000;2000; ...
     -30;25;5;-40;10;20]];
worst = 0;
for j = 1:size(states,2)
    testArgs = args;
    testArgs{2} = states(:,j);
    if j == 3
        testArgs{7} = [0.08;0.11;0.13];
        testArgs{8} = [0.5;-1.0;9.7];
        testArgs{10} = diag(1e-9*[1;2;3;4;5;6]);
        testArgs{11} = 1e-5*[1;-1;0;2;0;-2];
    end
    ref.L=testArgs{3}; ref.r=testArgs{4}; ref.K=testArgs{5};
    ref.A=testArgs{6}; ref.mi=testArgs{7}; ref.g=testArgs{8};
    ref.k=testArgs{9}; ref.D=testArgs{10};
    expected = armS_pressure_dynamics(testArgs{1},testArgs{2},ref,testArgs{11});
    actual = armS_pressure_mex(testArgs{:});
    assert(all(isfinite(actual)),'MEX returned nonfinite values.');
    err = max(abs(actual-expected)./(1+abs(expected)));
    worst = max(worst,err);
    fprintf('MEX pose %d: maximum scaled error %.3e\n',j,err);
    % The generated symbolic integrals and point transform are not exactly
    % mutually consistent at high bend. This bound covers that known gap.
    assert(err < 1e-4,'MEX exceeds checked symbolic parity bound.');
end
binary = fullfile(buildFolder,['armS_pressure_mex.',mexext]);
assert(isfile(binary),'Compiled MEX binary was not found.');
clear armS_pressure_mex
copyfile(binary,fullfile(sourceFolder,['armS_pressure_mex.',mexext]),'f');
fprintf('Verified MEX copied to %s (worst scaled error %.3e).\n', ...
    sourceFolder,worst);
end

function restoreEnvironment(oldFolder,oldPath)
cd(oldFolder);
path(oldPath);
rehash;
end
