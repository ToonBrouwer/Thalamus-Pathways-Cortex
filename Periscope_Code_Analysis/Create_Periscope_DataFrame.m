%% Create Delta F over F dataset
%add paths
addpath(genpath('C:\Users\motot\Documents\MATLAB\Periscope_Code_Analysis\Functions')) %helper functions
addpath(genpath('C:\Users\motot\Documents\MATLAB\double_periscope')) % Periscope videos


%% create dataset
% List of animals
animals = {'ane', 'awa'};  % Add as many as needed
num_recordings = 20;  % Each animal has 20 recordings
Fs = 100;  % Sampling rate in frames per second
trial_durations = [4, 4, 4];  % Duration of each trial (in seconds)
baseline_durations = [0.8, 0.8, 0.8];  % Baseline duration (in seconds)

% Initialize the structure to store results
all_data = struct();
masks = struct();  % Structure to hold masks for each animal

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
%% create dataframe
% Define the list of animals, trials, and flexible processing options
animals = {'ane', 'awa'};  % Add as many animals as needed
num_trials = 3;  % Number of trials
Fs = 100;  % Sampling rate in frames per second
trial_durations = [4, 4, 4];  % Trial durations in seconds
baseline_durations = [0.6, 0.6, 0.6];  % Baseline durations in seconds
num_recordings = 20;

% Parameters for flexible dataset creation
num_animals_to_process = 2;  % Number of animals to process
num_recordings_to_process = 20;  % Number of recordings per animal to process
include_roi1 = true;  % Option to include ROI1 in the dataset
include_roi2 = true;  % Option to include ROI2 in the dataset

% Initialize the structure to store results
all_data = struct();

% Loop over each animal (limited to 'num_animals_to_process')
for aIdx = 1:min(num_animals_to_process, length(animals))
    animal_name = animals{aIdx};
    
    % Initialize fields for the current animal
    all_data.(animal_name) = struct();
    
    % Loop over each recording (limited to 'num_recordings_to_process')
    for recIdx = 1:min(num_recordings_to_process, num_recordings)
        % Generate filename
        video_file = sprintf('%s_recording%d.avi', animal_name, recIdx);
        
        try
            % Load video and convert to grayscale
            frames = load_video_as_grayscale(video_file);
            
            % Compute ΔF/F for each trial within each ROI
            deltaF_F = compute_deltaF_F(frames, Fs, trial_durations, baseline_durations, masks.(animal_name));
            
            % Initialize structure fields for this recording and store signals
            for trialIdx = 1:num_trials
                % Create a structure for the current trial and recording
                trial_data = struct();
                
                % Optionally include ROI1 data
                if include_roi1
                    trial_data.ROI1 = deltaF_F{trialIdx, 1};
                end
                
                % Optionally include ROI2 data
                if include_roi2
                    trial_data.ROI2 = deltaF_F{trialIdx, 2};
                end
                
                % Store the trial data in the 'all_data' structure
                all_data.(animal_name).(['recording' num2str(recIdx)]).(['trial' num2str(trialIdx)]) = trial_data;
            end
            
            % Clear variables for the current recording
            clear frames deltaF_F;
            fprintf('Successfully processed %s \n. ', video_file);
            
        catch ME
            % If the file is missing or there's an error, print a warning and continue
            warning('Could not process %s. Error: %s', video_file, ME.message);
            
            % Initialize empty fields to keep the structure consistent
            for trialIdx = 1:num_trials
                all_data.(animal_name).(['recording' num2str(recIdx)]).(['trial' num2str(trialIdx)]) = struct( ...
                    'ROI1', [], ...
                    'ROI2', [] ...
                );
            end
        end
        
        % Clear variables that were set inside the loop but are no longer needed
        clear video_file;
    end
end

%% Save the results to a .mat file for future use
save('Combined_data_per_region.mat', 'all_data', 'masks');
