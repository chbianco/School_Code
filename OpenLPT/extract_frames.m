% EXTRACT_FRAMES  Extract all frames from a video file into a new folder.
%
% USAGE:
%   Run the script — a file picker will open.
%
% OUTPUT:
%   A folder named after the video file (without extension) is created in
%   the same directory as the video. Frames are saved as greyscale TIFs
%   named frame_00001.tif, frame_00002.tif, etc.

% ── Select video file ─────────────────────────────────────────────────────────
[fname, fpath] = uigetfile( ...
    {'*.mov;*.avi;*.mp4;*.mj2;*.mjpeg;*.m4v','Video Files'; ...
     '*.*','All Files'}, ...
    'Select a video file');

if isequal(fname, 0)
    disp('No file selected.');
    return
end

videoPath = fullfile(fpath, fname);
[~, videoName, ~] = fileparts(fname);

% ── Create output folder ──────────────────────────────────────────────────────
outFolder = fullfile(fpath, videoName);
if ~exist(outFolder, 'dir')
    mkdir(outFolder);
    fprintf('Created folder: %s\n', outFolder);
else
    fprintf('Folder already exists: %s\n', outFolder);
end

% ── Read and export frames ────────────────────────────────────────────────────
v = VideoReader(videoPath);
nFrames = floor(v.Duration * v.FrameRate);

fprintf('Video     : %s\n', fname);
fprintf('Duration  : %.2f s\n', v.Duration);
fprintf('Frame rate: %.2f fps\n', v.FrameRate);
fprintf('Frames    : ~%d\n\n', nFrames);

count = 0;
while hasFrame(v)
    frame = readFrame(v);

    % Take channel 1 — avoids any weighted blend artefacts.
    % For a true greyscale source all channels are identical anyway.
    if size(frame, 3) > 1
        imgGray = frame(:,:,1);
    else
        imgGray = frame;
    end

    count = count + 1;
    outPath = fullfile(outFolder, sprintf('frame_%05d.tif', count));
    imwrite(imgGray, outPath);

    if mod(count, 100) == 0
        fprintf('  Written %d / ~%d frames...\n', count, nFrames);
    end
end

fprintf('\nDone. %d frames saved to:\n  %s\n', count, outFolder);