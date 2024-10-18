function deltaF_F = compute_deltaF_F(frames, Fs, trial_durations, baseline_durations, masks)
% compute_deltaF_F
% This function computes the ΔF/F for each trial within specified ROIs
% while retaining the temporal dimension. Pixels outside the ROI mask
% are set to NaN, and ΔF/F is computed ignoring NaN values.
%
% Inputs:
% frames: matrix of grayscale frames (pixelheight x pixelwidth x nframes)
% Fs: sampling rate in frames per second (e.g., 100 Hz)
% trial_durations: array specifying the duration (in seconds) of each trial
% baseline_durations: array specifying the baseline period (in seconds) for each trial
% masks: cell array of binary masks specifying the ROIs (e.g., {mask1, mask2})
%
% Outputs:
% deltaF_F: cell array where each cell contains the ΔF/F time series for each ROI

numTrials = length(trial_durations);  % Number of trials
trial_durations_frames = trial_durations * Fs;  % Convert trial durations to frames
baseline_durations_frames = baseline_durations * Fs;  % Convert baseline durations to frames
numROIs = length(masks);  % Number of ROIs (e.g., 2 for Circle 1 and Circle 2)

% Pre-allocate a cell array to store ΔF/F time series for each trial and ROI
deltaF_F = cell(numTrials, numROIs);

% Convert frames to double for calculations
frames = double(frames);  % Ensure frames are double for computation

% Initialize the start frame for the first trial
trial_start = 1;

% Loop through each trial to compute ΔF/F
for trialIdx = 1:numTrials
    % Define the start and end of the current trial
    trial_end = trial_start + trial_durations_frames(trialIdx) - 1;

    % Extract trial frames
    trial_frames = frames(:, :, trial_start:trial_end);

    % Extract baseline frames for the trial
    baseline_frames = frames(:, :, trial_start:(trial_start + baseline_durations_frames(trialIdx) - 1));

    % Loop through each ROI and compute ΔF/F time series
    for roiIdx = 1:numROIs
        % Get the mask for the current ROI
        mask = masks{roiIdx};

        mask_3D = repmat(mask, [1, 1, size(baseline_frames, 3)]);  % Replicate the 2D mask across frames
        baseline_frames(~mask_3D) = NaN;  % Apply NaN outside the mask
        trial_frames(~mask_3D) = NaN;  % Apply NaN outside the mask


        % Compute the mean baseline (F0) within the ROI, ignoring NaNs
        F0 = mean(baseline_frames(mask_3D), 'omitnan');  % Mean over all valid pixels in ROI

        % Compute ΔF/F for the current trial across all frames
        deltaF_F_signal = (trial_frames - F0) ./ F0;  % Compute ΔF/F

        % Average ΔF/F across the pixels inside the mask for each frame (ignoring NaNs)
        deltaF_F{trialIdx, roiIdx} = squeeze(mean(deltaF_F_signal, [1, 2], 'omitnan'));
    end

    % Update trial_start for the next trial
    trial_start = trial_end + 1;
end
end
