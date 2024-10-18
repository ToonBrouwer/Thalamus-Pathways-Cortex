%% Analysis script over trials
addpath(genpath("C:\Users\motot\Documents\MATLAB\Periscope_Code_Analysis\Functions"))
 load("Combined_data_per_region.mat"); %load data

%% Create plot of Event related traces
% Define the list of animals and trials
animals = {'ane', 'awa'};  % Add as many animals as needed
num_trials = 3;  % Number of trials (e.g., 3 for control, flash, and tone)
Fs = 100;  % Sampling rate in frames per second
frame_range = 70:170;  % Focus on frames 70 to 170 (300 ms before and 300 ms after)
num_recordings = 20;

% Time vector for plotting: -300 ms to 700 ms
time_ms = (-300:10:700);  % 10 ms per frame, covering 100 frames

% Loop through each animal and plot the data
for aIdx = 1:length(animals)
    animal_name = animals{aIdx};
    
    % Create a new figure for this animal
    figure;
    
    % Loop over each trial type and create subplots
    for trialIdx = 1:num_trials
        % Initialize empty arrays for storing data across recordings
        trial_data_ROI1 = [];
        trial_data_ROI2 = [];
        
        % Check for missing trials
        has_data = false;  % Flag to track if data is found for this trial

        % Collect data across recordings
        for recIdx = 1:num_recordings
            try
                % Try to access the trial data for ROI1 and ROI2
                trial_ROI1 = all_data.(animal_name).(['recording' num2str(recIdx)]).(['trial' num2str(trialIdx)]).ROI1(frame_range);
                trial_ROI2 = all_data.(animal_name).(['recording' num2str(recIdx)]).(['trial' num2str(trialIdx)]).ROI2(frame_range);
                
                % If data is found, append it to the arrays
                trial_data_ROI1 = [trial_data_ROI1; trial_ROI1'];
                trial_data_ROI2 = [trial_data_ROI2; trial_ROI2'];
                has_data = true;  % Set flag to true when data is found
                
            catch
                % If trial is missing, skip this recording
                warning('Missing trial %d for recording %d of animal %s.', trialIdx, recIdx, animal_name);
                continue;
            end
        end
        
        % If no data is found for this trial, skip plotting
        if ~has_data
            warning('No data found for trial %d of animal %s. Skipping plot.', trialIdx, animal_name);
            continue;
        end
        
        % Compute the mean and SEM for both ROIs
        mean_ROI1 = mean(trial_data_ROI1, 1);
        mean_ROI2 = mean(trial_data_ROI2, 1);
        sem_ROI1 = std(trial_data_ROI1, 0, 1) / sqrt(size(trial_data_ROI1, 1));  % SEM for ROI1
        sem_ROI2 = std(trial_data_ROI2, 0, 1) / sqrt(size(trial_data_ROI2, 1));  % SEM for ROI2
        
        % Create subplot for this trial
        subplot(3, 1, trialIdx);
        
        % Plot mean and SEM for ROI1
        shadedErrorBar(time_ms, mean_ROI1, sem_ROI1, 'lineprops', {'r', 'DisplayName', 'Circle 1'});
        hold on;
        
        % Plot mean and SEM for ROI2
        shadedErrorBar(time_ms, mean_ROI2, sem_ROI2, 'lineprops', {'b', 'DisplayName', 'Circle 2'});
        
        % Add labels, title, and legend
        xlabel('Time (ms)');
        ylabel('ΔF/F Signal');
        title(sprintf('Animal: %s, Trial %d', animal_name, trialIdx));
        legend('Circle 1', 'Circle 2');
        
        % Add stimulus onset and offset lines for Trials 2 and 3
        xline(0, '--k', 'LineWidth', 1.5);  % Stimulus onset
        xline(400, '--k', 'LineWidth', 1.5);  % Stimulus offset
    end
    
end
