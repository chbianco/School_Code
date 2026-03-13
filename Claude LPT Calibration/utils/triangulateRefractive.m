function X = triangulateRefractive(observations, camPoses, Ks, distCoeffs, ...
                                    wallGeoms, n, cfg)
% TRIANGULATEREFRACTIVE  Reconstruct a 3-D point from multi-camera
%                        observations, accounting for refraction.
%
%   Inputs:
%     observations  [nCams×2]  image observations [u v] per camera
%                              (NaN rows for cameras that didn't see the point)
%     camPoses      {nCams}   struct array .R, .t
%     Ks            {nCams}   intrinsic matrices
%     distCoeffs    {nCams}   distortion coefficients
%     wallGeoms     {nCams}   wall geometry structs
%     n             struct    .air .glass .water
%     cfg           struct    .maxReprojErr
%
%   Output:
%     X   [3×1]  triangulated 3-D point (NaN if failed)
%                with attached field .reprojErrors [nCams×1]

    % Determine which cameras see this point
    visible = ~any(isnan(observations), 2);
    camIdx  = find(visible);
    nVis    = numel(camIdx);

    if nVis < 2
        X = struct('pos', [NaN;NaN;NaN], 'reprojErrors', NaN(size(observations,1),1));
        return;
    end

    % Initial estimate via standard DLT on undistorted + "straightened" rays
    % (ignoring refraction for initialisation only)
    X0 = dltInitialise(observations(camIdx,:), camPoses(camIdx), Ks(camIdx), distCoeffs(camIdx));

    % Nonlinear refinement minimising sum of squared reprojection errors
    cost = @(x) refractionReprojCost(x, observations, camIdx, ...
                                      camPoses, Ks, distCoeffs, wallGeoms, n);

    opts   = optimset('Display','off','TolX',1e-9,'TolFun',1e-9,'MaxIter',200);
    X_opt  = fminsearch(cost, X0, opts);

    % Compute per-camera reprojection errors
    nCams       = size(observations,1);
    reprojErrors = NaN(nCams,1);
    for c = camIdx'
        uv_proj = projectPointRefractive(X_opt, camPoses{c}, Ks{c}, ...
                                          distCoeffs{c}, wallGeoms{c}, n);
        uv_obs  = observations(c,:)';
        reprojErrors(c) = norm(uv_proj - uv_obs);
    end

    X = struct('pos', X_opt(:), 'reprojErrors', reprojErrors);
end

