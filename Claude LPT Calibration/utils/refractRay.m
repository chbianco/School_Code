function ray_out = refractRay(rayOrigin, rayDir, wallNormal, wallPoint, ...
                               glassThick, n_air, n_glass, n_water)

    rayDir     = rayDir(:)     / norm(rayDir);
    wallNormal = wallNormal(:) / norm(wallNormal);
    wallPoint  = wallPoint(:);
    rayOrigin  = rayOrigin(:);

    outerPoint = wallPoint + glassThick * wallNormal;
    
    
    if dot(wallNormal, rayDir) > 0
        wallNormal = -wallNormal;
    end

    [hitAir, tAir] = rayPlaneIntersect(rayOrigin, rayDir, wallNormal, outerPoint);
    if isnan(tAir) || tAir < 0
        ray_out = struct('origin',[],'dir',[],'valid',false,'hitAir',[],'hitGlass',[]);
        return;
    end

    dir_glass = snellRefract(rayDir, wallNormal, n_air, n_glass);
    if isempty(dir_glass)
        ray_out = struct('origin',[],'dir',[],'valid',false,'hitAir',hitAir,'hitGlass',[]);
        return;
    end

    [hitGlass, tGlass] = rayPlaneIntersect(hitAir, dir_glass, wallNormal, wallPoint);
    if isnan(tGlass) || tGlass < 0
        ray_out = struct('origin',[],'dir',[],'valid',false,'hitAir',hitAir,'hitGlass',[]);
        return;
    end

    dir_water = snellRefract(dir_glass, wallNormal, n_glass, n_water);
    if isempty(dir_water)
        ray_out = struct('origin',[],'dir',[],'valid',false,'hitAir',hitAir,'hitGlass',hitGlass);
        return;
    end

    ray_out = struct('origin',  hitGlass(:), ...
                     'dir',     dir_water(:), ...
                     'valid',   true, ...
                     'hitAir',  hitAir(:), ...
                     'hitGlass',hitGlass(:));
end