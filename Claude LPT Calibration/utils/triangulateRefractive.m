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

    % Project seed into water volume using wall geometry
    % Find the water-side extent of each wall and clip accordingly
    for wi = 1:numel(wallGeoms)
        nm = wallGeoms{wi}.normal(:);
        pt = wallGeoms{wi}.point(:);
        % Distance of X0 from this wall (positive = water side)
        d = dot(nm, X0(:) - pt);
        if d > 0
            % X0 is on air side of this wall — project it just inside
            X0 = X0(:) - (d + 0.01) * nm;
        end
    end
    X0 = X0(:)';

    % Nonlinear refinement minimising sum of squared reprojection errors
    cost = @(x) refractionReprojCost(x, observations, camIdx, ...
                                      camPoses, Ks, distCoeffs, wallGeoms, n);

    % Use fmincon to keep search within water volume
    lb_opt = zeros(1,3); ub_opt = zeros(1,3);
    for wi = 1:numel(wallGeoms)
        nm = wallGeoms{wi}.normal(:)';
        pt = wallGeoms{wi}.point(:)';
        for dim = 1:3
            if nm(dim) > 0.5
                ub_opt(dim) = pt(dim);
            elseif nm(dim) < -0.5
                lb_opt(dim) = pt(dim);
        end
    end
    end
    % Add margin
    lb_opt = lb_opt - 0.05;
    ub_opt = ub_opt + 0.05;
    
    opts = optimoptions('fmincon', 'Display', 'off', ...
        'TolX', 1e-9, 'TolFun', 1e-9, 'MaxIterations', 200);
    X_opt = fmincon(cost, X0, [], [], [], [], lb_opt, ub_opt, [], opts);
    X_opt = X_opt(:);  % force to column vector

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

function cost = refractionReprojCost(X, observations, camIdx, ...
                                      camPoses, Ks, distCoeffs, wallGeoms, n)
    cost = 0;
    for c = camIdx'
        uv_proj = projectPointRefractive(X(:), camPoses{c}, Ks{c}, ...
                                          distCoeffs{c}, wallGeoms{c}, n);
        uv_obs  = observations(c,:)';
        if any(isnan(uv_proj))
            % Penalty that grows with distance from volume centre
            cost = cost + 1e6 + 1e4 * norm(X);
        else
            cost = cost + sum((uv_proj - uv_obs).^2);
        end
    end
end

function X0 = dltInitialise(obs, poses, Ks, distCoeffs)
% DLT-based initialisation (no refraction) for triangulation seed.
    nVis = numel(poses);
    A    = zeros(2*nVis, 3);
    b    = zeros(2*nVis, 1);
    for i = 1:nVis
        K      = Ks{i};
        R      = poses{i}.R;
        t      = poses{i}.t;
        uv     = obs(i,:);
        % Undo distortion
        xn = (uv(1) - K(1,3)) / K(1,1);
        yn = (uv(2) - K(2,3)) / K(2,2);
        % Bearing vector in world coords
        bear = R' * [xn; yn; 1];
        bear = bear / norm(bear);
        % Camera centre
        C = -R' * t;
        % Two equations from cross product
        A(2*i-1, :) = [0, -bear(3), bear(2)];
        A(2*i,   :) = [bear(3), 0, -bear(1)];
        b(2*i-1)    = -dot([0,-bear(3),bear(2)], C);
        b(2*i)      = -dot([bear(3),0,-bear(1)], C);
    end
    X0 = A \ b;
end
