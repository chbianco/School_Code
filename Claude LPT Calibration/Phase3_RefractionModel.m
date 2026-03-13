% =========================================================================
% Phase3_RefractionModel.m  (utility — not run directly)
% =========================================================================
% This file defines the core refractive ray-tracing functions used by
% Phase4 (bundle adjustment) and Phase5 (reconstruction).
%
% PHYSICAL MODEL:
%   Each camera is in air.  Its optical axis passes through a flat glass
%   wall before entering the water volume.  The wall is described by a
%   unit outward normal and a point on the inner (water-side) face.
%
%   A ray from camera centre C in direction d_air must be refracted at:
%     Interface 1:  air   → glass   (n1=1.000, n2=1.520)
%     Interface 2:  glass → water   (n1=1.520, n2=1.333)
%
%   Snell's law (vector form):
%     n1*(d × n̂) = n2*(d' × n̂)
%     d' = (n1/n2)*d + (n1/n2*cos(θ_i) - cos(θ_t))*n̂
%     where cos(θ_i) = -d·n̂,  cos(θ_t) = sqrt(1 - (n1/n2)²*(1-cos²(θ_i)))
%
% FUNCTIONS (all can be called from other scripts):
%   ray  = refractRay(rayOrigin, rayDir, wallNormal, wallPoint, ...
%                     glassThick, n_air, n_glass, n_water)
%          Returns the final ray direction + origin INSIDE the water
%
%   uv   = projectPointRefractive(X_world, camPose, intrinsics, wallGeom, refractIdx)
%          Projects a 3-D world point through the full refractive model to
%          image coordinates.  Used in the bundle adjustment cost function.
%
%   X    = triangulateRefractive(observations, camPoses, intrinsics, wallGeoms, refractIdx)
%          Multi-camera triangulation respecting refraction.
%
%   dP   = refractionJacobian(X_world, camPose, intrinsics, wallGeom, refractIdx)
%          Numerical Jacobian of projection w.r.t. world point (for uncertainty).
% =========================================================================


% =========================================================================
function ray_out = refractRay(rayOrigin, rayDir, wallNormal, wallPoint, ...
                               glassThick, n_air, n_glass, n_water)
% REFRACTRAY  Trace a ray from air through a flat glass wall into water.
%
%   Inputs:
%     rayOrigin  [3×1]  ray start point in world coords (in air) [m]
%     rayDir     [3×1]  unit direction of ray in air
%     wallNormal [3×1]  outward unit normal of glass wall (pointing INTO air)
%     wallPoint  [3×1]  point on inner (water-side) face of glass [m]
%     glassThick  scalar  wall thickness [m]
%     n_air, n_glass, n_water  refractive indices
%
%   Output:
%     ray_out  struct with fields:
%       .origin    [3×1]  ray origin on water-side of glass
%       .dir       [3×1]  unit ray direction in water
%       .valid     logical  false if total internal reflection
%       .hitAir    [3×1]  intersection with outer (air-side) wall face
%       .hitGlass  [3×1]  intersection with inner (water-side) wall face

    rayDir    = rayDir    / norm(rayDir);
    wallNormal = wallNormal / norm(wallNormal);

    % Outer (air-side) face is offset from wallPoint by glassThick along normal
    outerPoint = wallPoint + glassThick * wallNormal;

    % ---- Interface 1: air → glass  (at outer face) ----------------------
    [hitAir, tAir] = rayPlaneIntersect(rayOrigin, rayDir, wallNormal, outerPoint);
    if tAir < 0
        % Ray is going away from wall — use -normal to detect approach from other side
        % This handles cameras that view through the floor or opposite wall
        [hitAir, tAir] = rayPlaneIntersect(rayOrigin, rayDir, -wallNormal, outerPoint);
        if tAir < 0
            ray_out = struct('origin',[],'dir',[],'valid',false,'hitAir',[],'hitGlass',[]);
            return;
        end
        wallNormal = -wallNormal;  % flip so refraction is computed correctly
    end

    dir_glass = snellRefract(rayDir, -wallNormal, n_air, n_glass);
    if isempty(dir_glass)
        ray_out = struct('origin',[],'dir',[],'valid',false,'hitAir',hitAir,'hitGlass',[]);
        return;
    end

    % ---- Interface 2: glass → water  (at inner face = wallPoint) ---------
    [hitGlass, tGlass] = rayPlaneIntersect(hitAir, dir_glass, wallNormal, wallPoint);
    if tGlass < 0
        ray_out = struct('origin',[],'dir',[],'valid',false,'hitAir',hitAir,'hitGlass',[]);
        return;
    end

    dir_water = snellRefract(dir_glass, -wallNormal, n_glass, n_water);
    if isempty(dir_water)
        ray_out = struct('origin',[],'dir',[],'valid',false,'hitAir',hitAir,'hitGlass',hitGlass);
        return;
    end

    ray_out = struct('origin', hitGlass, 'dir', dir_water, ...
                     'valid',  true, ...
                     'hitAir', hitAir, 'hitGlass', hitGlass);
end


% =========================================================================
function uv = projectPointRefractive(X_world, camPose, K, distCoeffs, wallGeom, n)
% PROJECTPOINTREFRACTIVE  Project a 3-D world point to image coords through
%                         the full refractive model.
%
%   Inputs:
%     X_world  [3×1]  point in world coordinates [m]
%     camPose  struct  .R [3×3], .t [3×1]  world-to-camera transform
%                      (such that X_cam = R*X_world + t)
%     K        [3×3]  camera intrinsic matrix
%     distCoeffs [1×5]  [k1 k2 p1 p2 k3]
%     wallGeom struct  .normal, .point, .thickness  (world coords)
%     n        struct  .air, .glass, .water
%
%   Output:
%     uv  [2×1]  image coordinates [u; v] in pixels (NaN if invalid)

    % Camera centre in world coords
    C = -camPose.R' * camPose.t;

    % Direction from camera to world point (in world coords)
    d_air = X_world - C;
    d_air = d_air / norm(d_air);

    % Refract through wall
    ray = refractRay(C, d_air, wallGeom.normal, wallGeom.point, ...
                     wallGeom.thickness, n.air, n.glass, n.water);

    if ~ray.valid
        uv = [NaN; NaN];
        return;
    end

    % Find where the in-water ray hits the plane through X_world
    % perpendicular to its propagation (or just use the refracted origin
    % and direction to project, accounting for the offset)
    %
    % We compute the "virtual" point: where the refracted ray in water
    % would arrive at X_world's depth, then project THAT point through
    % the standard pinhole model.

    % Express the water-side ray in camera coordinates
    origin_cam = camPose.R * ray.origin + camPose.t;
    dir_water_cam = camPose.R * ray.dir;

    % We need the intersection of the water ray with the plane
    % through X_world perpendicular to the camera Z axis.
    X_cam = camPose.R * X_world + camPose.t;

    % Parameterise: P(s) = origin_cam + s * dir_water_cam
    % Constraint: P_z(s) = X_cam(3)   (same depth)
    if abs(dir_water_cam(3)) < 1e-10
        uv = [NaN; NaN];
        return;
    end
    s      = (X_cam(3) - origin_cam(3)) / dir_water_cam(3);
    P_cam  = origin_cam + s * dir_water_cam;

    % Normalised image coords
    x_n    = P_cam(1) / P_cam(3);
    y_n    = P_cam(2) / P_cam(3);

    % Apply lens distortion
    [x_d, y_d] = applyDistortion(x_n, y_n, distCoeffs);

    % Apply intrinsic matrix
    u = K(1,1) * x_d + K(1,3);
    v = K(2,2) * y_d + K(2,3);

    uv = [u; v];
end


% =========================================================================
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


% =========================================================================
function J = refractionJacobian(X_world, camPose, K, distCoeffs, wallGeom, n)
% REFRACTIONJACOBIAN  Numerical 2×3 Jacobian d(uv)/d(X_world).
%                     Used for uncertainty propagation.

    eps  = 1e-6;
    uv0  = projectPointRefractive(X_world, camPose, K, distCoeffs, wallGeom, n);
    J    = zeros(2,3);
    for i = 1:3
        dX       = zeros(3,1); dX(i) = eps;
        uv_plus  = projectPointRefractive(X_world+dX, camPose, K, distCoeffs, wallGeom, n);
        uv_minus = projectPointRefractive(X_world-dX, camPose, K, distCoeffs, wallGeom, n);
        J(:,i)   = (uv_plus - uv_minus) / (2*eps);
    end
end


% =========================================================================
%  PRIVATE HELPERS
% =========================================================================

function d_out = snellRefract(d_in, n_hat, n1, n2)
% Vector form of Snell's law.
% d_in    unit incident direction
% n_hat   unit surface normal pointing INTO the incident medium
% Returns d_out (unit refracted direction), or [] for TIR.

    d_in  = d_in  / norm(d_in);
    n_hat = n_hat / norm(n_hat);
    ratio = n1 / n2;
    cosI  = -dot(n_hat, d_in);
    sin2T = ratio^2 * (1 - cosI^2);

    if sin2T > 1.0   % total internal reflection
        d_out = [];
        return;
    end

    cosT  = sqrt(1 - sin2T);
    d_out = ratio * d_in + (ratio * cosI - cosT) * n_hat;
    d_out = d_out / norm(d_out);
end


function [P, t] = rayPlaneIntersect(rayOrigin, rayDir, planeNormal, planePoint)
% Compute intersection of ray with plane.
% Returns intersection point P and parameter t (P = rayOrigin + t*rayDir).
% t = NaN if ray is parallel to plane.

    denom = dot(planeNormal, rayDir);
    if abs(denom) < 1e-12
        P = NaN(3,1); t = NaN; return;
    end
    t = dot(planeNormal, planePoint - rayOrigin) / denom;
    P = rayOrigin + t * rayDir;
end


function [xd, yd] = applyDistortion(xn, yn, distCoeffs)
% Apply Brown-Conrady distortion model.
% distCoeffs = [k1 k2 p1 p2 k3]

    if numel(distCoeffs) < 5, distCoeffs(end+1:5) = 0; end
    k1=distCoeffs(1); k2=distCoeffs(2);
    p1=distCoeffs(3); p2=distCoeffs(4); k3=distCoeffs(5);

    r2 = xn.^2 + yn.^2;
    r4 = r2.^2;
    r6 = r2.^3;

    radial = 1 + k1*r2 + k2*r4 + k3*r6;
    xd = xn .* radial + 2*p1.*xn.*yn     + p2*(r2 + 2*xn.^2);
    yd = yn .* radial + p1*(r2 + 2*yn.^2) + 2*p2.*xn.*yn;
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

        % Two equations from cross product (ray × direction = 0)
        A(2*i-1, :) = [0, -bear(3), bear(2)];
        A(2*i,   :) = [bear(3), 0, -bear(1)];
        b(2*i-1)    = -dot([0,-bear(3),bear(2)], C);
        b(2*i)      = -dot([bear(3),0,-bear(1)], C);
    end

    X0 = A \ b;
end


function cost = refractionReprojCost(X, obs, camIdx, poses, Ks, distCoeffs, wallGeoms, n)
% Sum of squared reprojection errors for triangulation optimisation.
    cost = 0;
    for c = camIdx'
        uv_proj = projectPointRefractive(X(:), poses{c}, Ks{c}, ...
                                          distCoeffs{c}, wallGeoms{c}, n);
        uv_obs  = obs(c,:)';
        if ~any(isnan(uv_proj))
            cost = cost + sum((uv_proj - uv_obs).^2);
        end
    end
end
