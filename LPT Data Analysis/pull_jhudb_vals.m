% token
authkey = 'edu.jhu.pha.turbulence.testing-201406';  
% dataset (channel flow re_tau ~ 1000)
dataset =  'channel';
% variable of interest
variable = 'velocity';
% temporal interpolation
temporal_method = 'pchip'; 
% spatial interpolation
spatial_method = 'lag4';
% spatial data capture
spatial_operator  = 'field';

time = 1.0;

nx = 16;
ny = 16;
nz = 16;
n_points = nx * ny * nz;

points = zeros(n_points,3);

gridX = linspace(4*pi, 4*pi + 1.6, nx); 
gridZ = linspace(0, 1.5*pi + 1.6, nz);
gridY = linspace(-1, -1 + 0.16, ny);

for i = 1:nx
    for j = 1:ny
        for k = 1:nz
            points((i - 1) * ny * nz + (j - 1) * nz + k, 1) = gridX(i);  
            points((i - 1) * ny * nz + (j - 1) * nz + k, 2) = gridY(j);
            points((i - 1) * ny * nz + (j - 1) * nz + k, 3) = gridZ(k);
        end
    end
end

velocity = getData(authkey, dataset, variable, time, temporal_method, spatial_method, spatial_operator, points);
