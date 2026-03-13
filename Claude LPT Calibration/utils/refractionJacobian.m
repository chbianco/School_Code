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

