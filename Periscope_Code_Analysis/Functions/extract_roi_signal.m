function avg_signals = extract_roi_signal(deltaF_F, createROIs, predefined_masks)
% extract_roi_signal
% This function extracts the average ΔF/F signal for ROIs across all frames.
% The function supports interactive ROI creation or using predefined masks.
%
% Inputs:
% deltaF_F: matrix of ΔF/F values (pixelheight x pixelwidth x nframes)
% createROIs: boolean flag (true to create ROIs interactively, false to use predefined masks)
% predefined_masks: (optional) cell array of predefined masks to use (if createROIs is false)
%
% Outputs:
% avg_signals: a matrix of size [nROIs x nFrames], where each row corresponds to the average signal for an ROI

    numFrames = size(deltaF_F, 3);  % Number of frames in ΔF/F data

    if createROIs
        % Create ROIs interactively using drawcircle or other methods
        figure;
        imshow(var(deltaF_F, 3), []);  % Display mean image across all frames for ROI selection
        title('Draw the first circle ROI');
        h1 = drawcircle();  % Interactive ROI drawing
        title('Draw the second circle ROI');
        h2 = drawcircle();  % Interactive ROI drawing

        % Create binary masks for each ROI
        masks{1} = createMask(h1);
        masks{2} = createMask(h2);
    else
        % If predefined masks are not provided, throw an error
        if nargin < 3 || isempty(predefined_masks)
            error('Predefined masks must be provided when createROIs is set to false.');
        end
        masks = predefined_masks;
    end

    numROIs = length(masks);  % Number of ROIs

    % Pre-allocate the output matrix to store average signals for each ROI
    avg_signals = zeros(numROIs, numFrames);

    % Use matrix operations to extract the average signal in each ROI
    for roiIdx = 1:numROIs
        % Apply the mask across all frames at once and calculate the average signal
        masked_frames = deltaF_F .* masks{roiIdx};  % Multiply mask with frames
        total_pixels = sum(masks{roiIdx}(:));  % Total number of pixels in the mask
        avg_signals(roiIdx, :) = squeeze(sum(sum(masked_frames, 1), 2)) / total_pixels;
    end
end
