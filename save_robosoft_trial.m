function files = save_robosoft_trial(trialDir, commands, pose, metadata)
%SAVE_ROBOSOFT_TRIAL Save one three-PMA command/pose experiment.
% commands: table with t_s,p_cmd1_Pa,p_cmd2_Pa,p_cmd3_Pa.
% pose: table with t_s,tip_x_m,tip_y_m,tip_z_m,theta_rad,phi_rad.
% Times are seconds on a shared clock. Actual chamber pressure is not logged.
% Existing trial files are never overwritten.

arguments
    trialDir (1,1) string
    commands table
    pose table
    metadata (1,1) struct
end

commandNames = ["t_s", "p_cmd1_Pa", "p_cmd2_Pa", "p_cmd3_Pa"];
poseNames = ["t_s", "tip_x_m", "tip_y_m", "tip_z_m", ...
             "theta_rad", "phi_rad"];
validateLog(commands, commandNames, false, "commands");
validateLog(pose, poseNames, true, "pose");

requiredMeta = ["trial_id", "trial_role", "command_unit", ...
                "pose_frame", "theta_definition", "phi_definition"];
missingMeta = requiredMeta(~isfield(metadata, cellstr(requiredMeta)));
if ~isempty(missingMeta)
    error("save_robosoft_trial:MissingMetadata", ...
        "Missing metadata: %s", strjoin(missingMeta, ", "));
end
if ~strcmpi(string(metadata.command_unit), "Pa")
    error("save_robosoft_trial:PressureUnit", ...
        "Convert commanded pressures to Pa before export.");
end

trialDir = char(trialDir);
files = struct("commands", fullfile(trialDir, "commands.csv"), ...
               "pose", fullfile(trialDir, "pose.csv"), ...
               "metadata", fullfile(trialDir, "trial.json"));
paths = struct2cell(files);
if any(cellfun(@(p) isfile(p) || isfolder(p), paths))
    error("save_robosoft_trial:Exists", ...
        "A trial output already exists; choose a new trial directory.");
end
if ~isfolder(trialDir)
    mkdir(trialDir);
end

metadata.schema_version = 1;
metadata.actual_chamber_pressure_measured = false;
metadata.command_samples = height(commands);
metadata.pose_samples = height(pose);

writetable(commands, files.commands);
writetable(pose, files.pose);
fid = fopen(files.metadata, "w");
if fid < 0
    error("save_robosoft_trial:WriteFailed", ...
        "Could not create %s", files.metadata);
end
cleanup = onCleanup(@() fclose(fid));
fwrite(fid, jsonencode(metadata), "char");
end

function validateLog(logTable, requiredNames, allowNaN, label)
names = string(logTable.Properties.VariableNames);
missing = requiredNames(~ismember(requiredNames, names));
if ~isempty(missing)
    error("save_robosoft_trial:MissingColumns", ...
        "%s lacks columns: %s", label, strjoin(missing, ", "));
end
if height(logTable) < 2
    error("save_robosoft_trial:TooShort", ...
        "%s must have at least two rows.", label);
end
t = logTable.t_s;
if ~isnumeric(t) || ~isvector(t) || any(~isfinite(t)) || ...
        any(diff(t) <= 0)
    error("save_robosoft_trial:Time", ...
        "%s.t_s must be finite and strictly increasing.", label);
end
for name = requiredNames(2:end)
    value = logTable.(name);
    if ~isnumeric(value) || ~isvector(value) || numel(value) ~= height(logTable) ...
            || ~isreal(value) || any(isinf(value)) ...
            || (~allowNaN && any(isnan(value)))
        error("save_robosoft_trial:Data", ...
            "%s.%s must be a numeric column; commands must be finite.", ...
            label, name);
    end
end
end
