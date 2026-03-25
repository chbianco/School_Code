  RADIUS = 10.0;
  METRIC = 0.5;
  NUM_CAMS = 4;
  OUTPUT_FILE = 'wand_points.csv';

  % --- Step 1: Find frames valid in ALL cameras ---
  common_frames = [];
  for cam = 1:NUM_CAMS
      d = detections{cam, 1};
      cam_valid_frames = [d([d.valid] == 1).frame];
      if cam == 1
          common_frames = cam_valid_frames;
      else
          common_frames = intersect(common_frames, cam_valid_frames);
      end
  end
  common_frames = sort(common_frames);
  fprintf('Found %d frames valid across all cameras.\n', numel(common_frames));

  % --- Step 2: Write CSV ---
  fid = fopen(OUTPUT_FILE, 'w');
  fprintf(fid, 'Frame,Camera,Status,PointIdx,X,Y,Radius,Metric\n');

  for fi = 1:numel(common_frames)
      frame_id = common_frames(fi);

      for cam = 1:NUM_CAMS
          d = detections{cam, 1};
          idx = find([d.frame] == frame_id, 1);

          pts = d(idx).pts;   % 2x2: [x1,y1; x2,y2]
          x_large = pts(1,1);  y_large = pts(1,2);  % row 1 = Large
          x_small = pts(2,1);  y_small = pts(2,2);  % row 2 = Small

          cam_id = cam - 1;   % 0-indexed

          % Small first (insertion order matters), then Large
          fprintf(fid, '%d,%d,Filtered_Small,0,%.6f,%.6f,%.1f,%.1f\n', ...
              frame_id, cam_id, x_small, y_small, RADIUS, METRIC);
          fprintf(fid, '%d,%d,Filtered_Large,1,%.6f,%.6f,%.1f,%.1f\n', ...
              frame_id, cam_id, x_large, y_large, RADIUS, METRIC);
      end
  end

  fclose(fid);
  fprintf('Wrote %d frames x %d cameras to %s\n', numel(common_frames), NUM_CAMS, OUTPUT_FILE);