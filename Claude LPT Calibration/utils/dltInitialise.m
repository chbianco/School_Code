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

