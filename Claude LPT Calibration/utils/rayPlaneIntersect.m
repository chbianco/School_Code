function [P, t] = rayPlaneIntersect(rayOrigin, rayDir, planeNormal, planePoint)
% Compute intersection of ray with plane.
% Returns intersection point P and parameter t (P = rayOrigin + t*rayDir).
% t = NaN if ray is parallel to plane.

    denom = dot(planeNormal, rayDir);
    if abs(denom) < 1e-12
        P = NaN(3,1); t = NaN; return;
    end
    t = dot(planeNormal, planePoint - rayOrigin) / denom;
    P = rayOrigin(:) + t * rayDir(:);
end

