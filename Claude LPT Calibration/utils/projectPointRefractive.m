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
    ray = refractRay(C(:), d_air(:), wallGeom.normal(:), wallGeom.point(:), ...
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
