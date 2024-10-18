function frames = load_video_as_grayscale(video_file)
%load_video_as_greyscale
% This function loads in a video at a given path as a matrix of grey-scaled
% frames.
% Inputs:
% video_file: path of the vidoe to read or filename if in CD 
% Outputs:
% frames: pixelheight x pixelwidth x nframes matrix of uint8 greyscale
% values

    % Load video data
    video_obj = VideoReader(video_file);

    % Pre-allocate grayscale frames structure
    numFrames = floor(video_obj.Duration * video_obj.FrameRate);  % Calculate number of frames
    frames = zeros(video_obj.Height, video_obj.Width, numFrames, 'uint8');  % Pre-allocate for grayscale

    % Load and convert frames to grayscale
    k = 1;
    while hasFrame(video_obj)
        rgbFrame = readFrame(video_obj);
        frames(:, :, k) = rgb2gray(rgbFrame);  % Convert the frame to grayscale and store
        k = k + 1;
    end
end
