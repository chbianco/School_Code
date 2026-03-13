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
