% Define the directory containing the .dat files
dataDir = 'data'; % Directory containing the .dat files
outputPatternT = 'T_output_*.dat'; % Pattern to match T .dat files
outputPatternS = 'S_output_*.dat'; % Pattern to match S .dat files

% Get a list of all .dat files for T and S
fileListT = dir(fullfile(dataDir, outputPatternT));
fileListS = dir(fullfile(dataDir, outputPatternS));
numFiles = length(fileListT);

% Ensure there are files to plot and the number of T and S files match
if numFiles == 0
    error('No .dat files found in the specified directory.');
elseif numFiles ~= length(fileListS)
    error('Mismatch between T and S files.');
end

% Create a video writer object for MP4 format
videoFileName = 'random_width16.mp4'; % Name of the output video file
video = VideoWriter(videoFileName, 'MPEG-4'); % Use MPEG-4 profile for MP4 format
video.FrameRate = 10; % Set the frame rate for the video
open(video); % Open the video file for writing

% Create a figure for the animation
figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]); % Fullscreen figure for larger animation
colormap('jet'); % Set a colormap for better visualization

% Loop over each file and plot the data for both T and S
for k = 1:numFiles
    % Construct the full file paths
    filePathT = fullfile(dataDir, fileListT(k).name);
    filePathS = fullfile(dataDir, fileListS(k).name);
    
    % Read the data from the .dat files for T and S
    dataT = load(filePathT);
    dataS = load(filePathS);
    
    % Automatically check and display the size of the data
    [numRowsT, numColsT] = size(dataT);
    [numRowsS, numColsS] = size(dataS);
    fprintf('Reading file %s with size %d x %d\n', fileListT(k).name, numRowsT, numColsT);
    fprintf('Reading file %s with size %d x %d\n', fileListS(k).name, numRowsS, numColsS);
    
    % Transpose the data to swap the X and Y axes
    dataT = dataT';
    dataS = dataS';
    
    % Plot the filled contour for T
    subplot(2, 1, 1); % First subplot (filled contour plot for T)
    contourf(dataT, 20, 'LineColor', 'none'); % Plot filled contour with 20 levels and no lines
    title(sprintf('Filled Contour Plot of T (Swapped Axes)', k));
    xlabel('Y'); % Y-axis is now horizontal
    ylabel('X'); % X-axis is now vertical
    colorbar; % Add a color bar for reference
    axis equal tight; % Adjust axis properties for better visualization
    
    % Plot the filled contour for S
    subplot(2, 1, 2); % Second subplot (filled contour plot for S)
    contourf(dataS, 20, 'LineColor', 'none'); % Plot filled contour with 20 levels and no lines
    title(sprintf('Filled Contour Plot of S (Swapped Axes)', k));
    xlabel('Y'); % Y-axis is now horizontal
    ylabel('X'); % X-axis is now vertical
    colorbar; % Add a color bar for reference
    axis equal tight; % Adjust axis properties for better visualization
    
    % Capture the current frame and write it to the video
    frame = getframe(gcf);
    writeVideo(video, frame);
    
    % Update the plot
    drawnow;
    
    % Pause for a brief moment to create the animation effect
    pause(0.01); % Adjust the pause duration as needed
end

% Close the video file
close(video);

disp('Animation completed and video saved as MP4.');
