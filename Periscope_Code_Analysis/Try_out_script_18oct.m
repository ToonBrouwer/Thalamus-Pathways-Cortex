% Define variables
animals = {'ane' 'awa'};  % Add as many animals as needed
num_trials = 3;  % Number of trials
Fs = 100;  % Sampling rate in frames per second
trial_durations = [4, 4, 4];  % Trial durations in seconds
baseline_durations = [0.5, 0.5, 0.5];  % Baseline durations in seconds
num_recordings = 10;  % Loop over recordings from 1 to 20

% Initialize the structure to store results across all recordings
all_data = struct();
%% Loop over each animal to define the masks
for aIdx = 1:length(animals)
    animal_name = animals{aIdx};
    
    % Load the first recording to help create the mask
    video_file = sprintf('%s_recording1.avi', animal_name);
    frames = load_video_as_grayscale(video_file);
    mean_image = mean(frames, 3);  % Compute the mean image across all frames to visualize

    % Show the image and let the user draw masks for this animal
    figure;
    imshow(mean_image, []);
    title(sprintf('Draw the first ROI for animal: %s', animal_name));
    h1 = drawcircle();  % Interactive ROI for first region
    title(sprintf('Draw the second ROI for animal: %s', animal_name));
    h2 = drawcircle();  % Interactive ROI for second region
    
    % Create binary masks for each ROI
    masks.(animal_name){1} = createMask(h1);  % Mask for first region
    masks.(animal_name){2} = createMask(h2);  % Mask for second region
    
    % Close the figure after mask selection
    close;
    clear frames mean_image h1 h2
end
%% Check fluorescence data
% Load the masks

% Preallocate the structure to store results across all recordings
all_data = struct();

for aIdx = 1:length(animals)
    animal_name = animals{aIdx};
      
    % Loop over each recording
    for recIdx = 1:num_recordings
        % Try to load video and process it
        video_file = sprintf('%s_recording%d.avi', animal_name, recIdx);
        
        try
            frames = load_video_as_grayscale(video_file);
            frames = double(frames);  % Convert to double for calculations

            % Initialize masks for this animal
            mask = masks.(animal_name);

            % Compute ΔF/F and raw fluorescence for this recording
            [deltaF_F, raw_fluorescence] = compute_deltaF_F_per_pixel(frames, Fs, trial_durations, baseline_durations, mask);

            % Store raw fluorescence and ΔF/F for this recording
            for trialIdx = 1:num_trials
                all_data.(animal_name).(['recording' num2str(recIdx)]).(['trial' num2str(trialIdx)]) = struct( ...
                    'ROI1_raw', raw_fluorescence{trialIdx, 1}, ...
                    'ROI2_raw', raw_fluorescence{trialIdx, 2}, ...
                    'ROI1_deltaF_F', deltaF_F{trialIdx, 1}, ...
                    'ROI2_deltaF_F', deltaF_F{trialIdx, 2} ...
                );
            end

            % Clear variables for the current recording
            clear frames deltaF_F raw_fluorescence;
            fprintf('Successfully processed %s\n', video_file);

        catch ME
            % Display warning and continue to the next recording
            warning('Could not process %s. Error: %s', video_file, ME.message);
            continue;
        end
    end
end

% Save the processed data
save('All_Combined_data_perpix.mat', 'all_data');


%% plot
% Load the data structure (previously saved)
 % load('All_Combined_data.mat');
% Define the animal you want to plot

animals = {'ane' 'awa'};  % Add as many animals as needed
num_trials = 3;  % Number of trials (e.g., 3 for control, flash, and tone)
Fs = 100;  % Sampling rate in frames per second
frame_range = 70:170;  % Focus on frames 70 to 170 (300 ms before and 300 ms after)
num_recordings = 10;  % Assume 8 recordings for now

% Time vector for plotting: -300 ms to 700 ms
time_ms = (-300:10:700);  % 10 ms per frame, covering 100 frames

% Loop through each animal
for aIdx = 1:length(animals)
    animal_name = animals{aIdx};
    
    % Plot raw fluorescence signals
    figure('Name', sprintf('%s - Raw Fluorescence', animal_name));
    
    % Loop over each trial type and create subplots for raw fluorescence
    for trialIdx = 1:num_trials
        subplot(3, 1, trialIdx);  % Subplot for each trial
        
        hold on;  % Allow multiple traces in one plot
        
        % Loop over each recording
        for recIdx = 1:num_recordings
            try
                % Extract raw fluorescence data for ROI1 and ROI2
                trial_ROI1 = all_data.(animal_name).(['recording' num2str(recIdx)]).(['trial' num2str(trialIdx)]).ROI1_raw(frame_range);
                trial_ROI2 = all_data.(animal_name).(['recording' num2str(recIdx)]).(['trial' num2str(trialIdx)]).ROI2_raw(frame_range);
                
                % Plot raw fluorescence for both ROIs
                hold on
                plot(time_ms, trial_ROI1, 'r');  % ROI 1
                % plot(time_ms, trial_ROI2, 'b');  % ROI 2
                
            catch
                % Handle missing trials by skipping the plot
                warning('Missing trial %d for recording %d of animal %s.', trialIdx, recIdx, animal_name);
                continue;
            end
        end
        
        % Add labels, title, and grid
        xlabel('Time (ms)');
        ylabel('Raw Fluorescence');
        title(sprintf('Animal: %s, Trial %d', animal_name, trialIdx));
        xline(0, '--k', 'LineWidth', 1.5);  % Stimulus onset
        xline(400, '--k', 'LineWidth', 1.5);  % Stimulus offset (for Trials 2 and 3)
        grid on;
    end
    
    %%Plot ΔF/F signals
    figure('Name', sprintf('%s - ΔF/F Signals', animal_name));
    
    % Loop over each trial type and create subplots for ΔF/F signals
    for trialIdx = 1:num_trials
        subplot(3, 1, trialIdx);  % Subplot for each trial
        
        hold on;  % Allow multiple traces in one plot
        
        % Loop over each recording
        for recIdx = 1:num_recordings
            try
                % Extract ΔF/F data for ROI1 and ROI2
                trial_ROI1 = all_data.(animal_name).(['recording' num2str(recIdx)]).(['trial' num2str(trialIdx)]).ROI1_deltaF_F(frame_range);
                trial_ROI2 = all_data.(animal_name).(['recording' num2str(recIdx)]).(['trial' num2str(trialIdx)]).ROI2_deltaF_F(frame_range);
                
                % Plot ΔF/F for both ROIs
                plot(time_ms, trial_ROI1, 'r');  % ROI 1
                plot(time_ms, trial_ROI2, 'b');  % ROI 2
                
            catch
                % Handle missing trials by skipping the plot
                warning('Missing trial %d for recording %d of animal %s.', trialIdx, recIdx, animal_name);
                continue;
            end
        end
        
        % Add labels, title, and grid
        xlabel('Time (ms)');
        ylabel('ΔF/F Signal');
        title(sprintf('Animal: %s, Trial %d', animal_name, trialIdx));
        xline(0, '--k', 'LineWidth', 1.5);  % Stimulus onset
        xline(400, '--k', 'LineWidth', 1.5);  % Stimulus offset (for Trials 2 and 3)
    end
end


function [deltaF_F, raw_fluorescence] = compute_deltaF_F(frames, Fs, trial_durations, baseline_durations, masks)
    % compute_deltaF_F
    % This function computes ΔF/F and raw fluorescence for multiple ROIs
    % Inputs:
    % frames: 3D matrix of video frames (height x width x time)
    % Fs: Sampling frequency (Hz)
    % trial_durations: array of trial durations (in seconds)
    % baseline_durations: array of baseline durations (in seconds)
    % masks: cell array of binary masks for each ROI
    %
    % Outputs:
    % deltaF_F: cell array of ΔF/F values per ROI
    % raw_fluorescence: cell array of raw fluorescence values per ROI

    numTrials = length(trial_durations);  % Number of trials
    trial_durations_frames = trial_durations * Fs;  % Convert trial durations to frames
    baseline_durations_frames = baseline_durations * Fs;  % Convert baseline durations to frames
    numROIs = length(masks);  % Number of ROIs (e.g., 2 for Circle 1 and Circle 2)

    % Pre-allocate cell arrays to store ΔF/F and raw fluorescence
    deltaF_F = cell(numTrials, numROIs);
    raw_fluorescence = cell(numTrials, numROIs);

    % Loop through each trial
    trial_start = 1;
    for trialIdx = 1:numTrials
        % Define the start and end of the current trial
        trial_end = trial_start + trial_durations_frames(trialIdx) - 1;

        % Extract trial frames
        trial_frames = frames(:, :, trial_start:trial_end);

        % Extract baseline frames
        baseline_frames = frames(:, :, trial_start:(trial_start + baseline_durations_frames(trialIdx) - 1));

        % Loop through each ROI
        for roiIdx = 1:numROIs
            % Get the mask for the current ROI
            mask = masks{roiIdx};

            % Apply the NaN approach to mask out irrelevant pixels
            masked_baseline_frames = baseline_frames;
            masked_trial_frames = trial_frames;
            masked_baseline_frames(~mask) = NaN;  % Mask out pixels outside ROI
            masked_trial_frames(~mask) = NaN;

            % Compute mean fluorescence signal across the pixels in the ROI (ignoring NaNs)
            mean_fluorescence = squeeze(mean(masked_trial_frames, [1, 2], 'omitnan'));
            mean_baseline_fluorescence = squeeze(mean(masked_baseline_frames, [1, 2], 'omitnan'));

            % Compute the baseline fluorescence (F0) by averaging baseline period (ignoring NaNs)
            F0 = mean(mean_baseline_fluorescence, 'omitnan');

            % Store raw fluorescence signal (mean pixel intensity per frame in ROI)
            raw_fluorescence{trialIdx, roiIdx} = mean_fluorescence;

            % Compute ΔF/F for the mean fluorescence signal in the ROI
            deltaF_F_signal = (mean_fluorescence - F0) / F0;

            % Store ΔF/F signal for this ROI
            deltaF_F{trialIdx, roiIdx} = deltaF_F_signal;
        end

        % Update trial_start for the next trial
        trial_start = trial_end + 1;
    end
end
function [deltaF_F, raw_fluorescence] = compute_deltaF_F_per_pixel(frames, Fs, trial_durations, baseline_durations, masks)
    % compute_deltaF_F_per_pixel
    % This function computes ΔF/F and raw fluorescence by first calculating ΔF/F per pixel and then averaging across pixels.
    % Inputs:
    % frames: 3D matrix of video frames (height x width x time)
    % Fs: Sampling frequency (Hz)
    % trial_durations: array of trial durations (in seconds)
    % baseline_durations: array of baseline durations (in seconds)
    % masks: cell array of binary masks for each ROI
    %
    % Outputs:
    % deltaF_F: cell array of ΔF/F values per ROI (computed by first calculating per pixel and then averaging)
    % raw_fluorescence: cell array of raw fluorescence values per ROI

    numTrials = length(trial_durations);  % Number of trials
    trial_durations_frames = trial_durations * Fs;  % Convert trial durations to frames
    baseline_durations_frames = baseline_durations * Fs;  % Convert baseline durations to frames
    numROIs = length(masks);  % Number of ROIs (e.g., 2 for Circle 1 and Circle 2)

    % Pre-allocate cell arrays to store ΔF/F and raw fluorescence
    deltaF_F = cell(numTrials, numROIs);
    raw_fluorescence = cell(numTrials, numROIs);

    % Loop through each trial
    trial_start = 1;
    for trialIdx = 1:numTrials
        % Define the start and end of the current trial
        trial_end = trial_start + trial_durations_frames(trialIdx) - 1;

        % Extract trial frames (height x width x time)
        trial_frames = frames(:, :, trial_start:trial_end);

        % Extract baseline frames (height x width x time)
        baseline_frames = frames(:, :, trial_start:(trial_start + baseline_durations_frames(trialIdx) - 1));

        % Loop through each ROI
        for roiIdx = 1:numROIs
            % Get the mask for the current ROI
            mask = masks{roiIdx};

            % Apply the NaN approach to mask out irrelevant pixels
            masked_baseline_frames = baseline_frames;
            masked_trial_frames = trial_frames;
            masked_baseline_frames(~mask) = NaN;  % Mask out pixels outside ROI
            masked_trial_frames(~mask) = NaN;

            % Compute baseline fluorescence (F0) per pixel
            F0_per_pixel = mean(masked_baseline_frames, 3, 'omitnan');  % Compute mean across time for baseline period

            % Compute ΔF/F for each pixel in the ROI across the trial duration
            deltaF_F_per_pixel = (masked_trial_frames - F0_per_pixel) ./ F0_per_pixel;

            % Compute the average ΔF/F across all pixels within the ROI at each time point
            deltaF_F_signal = squeeze(mean(deltaF_F_per_pixel, [1, 2], 'omitnan'));

            % Compute the mean raw fluorescence signal at each time point (spatial average within the ROI)
            raw_fluorescence_signal = squeeze(mean(masked_trial_frames, [1, 2], 'omitnan'));

            % Store raw fluorescence and ΔF/F for this trial and ROI
            raw_fluorescence{trialIdx, roiIdx} = raw_fluorescence_signal;
            deltaF_F{trialIdx, roiIdx} = deltaF_F_signal;
        end

        % Update trial_start for the next trial
        trial_start = trial_end + 1;
    end
end
