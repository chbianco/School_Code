function [u,v,w,p,x,y,z,nx,ny,nz  ] = load_channel_data( file_name )

tmp = load(file_name);

% number of grid points in x,y,z directions
nx = 100+1;
nz = 100+1;
ny = 193+1;

% equidistant grid spacing in x and z directions
x = linspace(0, 2*pi, nx);  x=x';
z = linspace(0, 1*pi, nz);  z=z';

% sort unstructured data into structured one using xyz coordinate info. 
% u in the end contains only the lower half of the channel.
% size(u) = [nx, ny/2, nz].

uuu = [tmp(:,1) tmp(:,2) tmp(:,3) tmp(:,4) ];
yyy = sortrows(uuu, [3 1 2]);
uuu = sortrows(uuu, [3 2 1]);
uu = uuu(:,4);
uu = reshape(uu, nx, ny, []);
u = uu(:, 1:ny/2, :);
u(:,1,:)=0;

% y coordinates: from wall to the close to the channel center line
y = yyy(1:ny/2, 2)+1;


vvv = [tmp(:,1) tmp(:,2) tmp(:,3) tmp(:,5) ];
vvv = sortrows(vvv, [3 2 1]);
vv = vvv(:,4);
vv = reshape(vv, nx, ny, []);
v = vv(:, 1:ny/2, :);
v(:,1,:)=0;

www = [tmp(:,1) tmp(:,2) tmp(:,3) tmp(:,6) ];
www = sortrows(www, [3 2 1]);
ww = www(:,4);
ww = reshape(ww, nx, ny, []);
w = ww(:, 1:ny/2, :);
w(:,1,:)=0;

ppp = [tmp(:,1) tmp(:,2) tmp(:,3) tmp(:,7) ];
ppp = sortrows(ppp, [3 2 1]);
pp = ppp(:,4);
pp = reshape(pp, nx, ny, []);
p = pp(:, 1:ny/2, :);


end

