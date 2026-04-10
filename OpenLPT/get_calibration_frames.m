% Script to extract first frame from all .mov files in a folder

% --- Select folder via dialog ---
folder_path = uigetdir(pwd, 'Select folder containing .mov files');

if folder_path == 0
    disp('No folder selected. Exiting.');
    return;
end

% --- Get all .mov files ---
mov_files = dir(fullfile(folder_path, '*.mov'));

if isempty(mov_files)
    fprintf('No .mov files found in: %s\n', folder_path);
    return;
end

fprintf('Found %d .mov file(s). Processing...\n', numel(mov_files));

% --- Process each file ---
for i = 1:numel(mov_files)
    video_name = mov_files(i).name;
    video_path = fullfile(folder_path, video_name);
    
    % Get name without extension
    [~, name_no_ext, ~] = fileparts(video_name);
    
    % Output filename
    output_filename = fullfile(folder_path, [name_no_ext '_frame.tif']);
    
    try
        % Read video and extract first frame
        v = VideoReader(video_path);
        first_frame = readFrame(v);
        
        % Save as .tif
        imwrite(first_frame, output_filename);
        fprintf('  [OK] %s -> %s\n', video_name, [name_no_ext '_frame.tif']);
        
    catch ME
        fprintf('  [ERROR] %s: %s\n', video_name, ME.message);
    end
end

fprintf('Done.\n');