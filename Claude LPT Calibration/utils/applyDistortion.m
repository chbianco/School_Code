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
