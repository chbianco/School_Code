function [Vol2U Vol2V Vol2S Vol2A Erange] = Xwire_compute_calibration(E1, E2, Vel, Theta)

% Calculate calibration surface for xwire using polnomial surface
%
% Two calibration modes - one in terms of speed and angle, the other in
% terms of U and V
%
% input is the raw calibration data from xwire_calibrate_acquire.m
%
% output is Vol2U and Vol2V - convert to U,V
%           Vol2S and Vol2A - convert to Speed and Angle
%
% Kenny Breuer, April 2025

PLOTTING = 1;

[Nangles Nspeed_cal] = size(E1);

%reduce the rank for data fitting
E1s = reshape(E1(1,:,:),1,Nangles*Nspeed_cal);
E2s = reshape(E2(1,:,:),1,Nangles*Nspeed_cal);
Vel = reshape(Vel,  1,Nangles*Nspeed_cal);
Ang = reshape(Theta,1,Nangles*Nspeed_cal);

fprintf('Caiibration data:\n')
fprintf('Velocity: %5.1f - %5.1f; Range: %5.2f [m/s]\n', min(Vel), max(Vel), max(Vel) - min(Vel));
fprintf('   Angle: %5.1f - %5.1f; Range: %5.2f [deg]\n', min(Ang), max(Ang), max(Ang) - min(Ang));
fprintf('      E1: %5.2f - %5.2f; Range: %5.2f [V]\n', min(E1s), max(E1s), max(E1s)-min(E1s));
fprintf('      E2: %5.2f - %5.2f; Range: %5.2f [V]\n', min(E2s), max(E2s), max(E2s)-min(E2s));

% Range of voltages
E1min = min(E1s);
E1max = max(E1s);
E2min = min(E2s);
E2max = max(E2s);
Erange = [E1min E1max E2min E2max];

% Calculate U and V of calibration velocities
Uvel = Vel.*cosd(Ang);
Vvel = Vel.*sind(Ang);

Angmin = min(Ang);
Angmax = max(Ang);

Vmin = min(Vvel);
Vmax = max(Vvel);

% calculate calibration library based on a 5th order polynomial
        Library = [ ...
            ones(1,length(E1s)); ...
            E1s; E2s; ...
            E1s.^2; E2s.^2; E1s.*E2s; ... 
            E1s.^3; E2s.^3; E1s.*E2s.^2; E1s.^2.*E2s; ...
            E1s.^4; E2s.^4; E1s.*E2s.^3; E1s.^3.*E2s; E1s.^2.*E2s.^2; ...
            E1s.^5; E2s.^5; E1s.*E2s.^4; E1s.^4.*E2s; E1s.^2.*E2s.^3; E1s.^3.*E2s.^2; ...
            ];

Library = Library';

A1 = Library\Vel';
Err_vel = std(100*(Library*A1 - Vel')./Vel');
A2 = Library\Ang';
Err_ang = std(Library*A2 - Ang');
A3 = Library\Uvel';
Err_Uvel = std(100*(Library*A3 - Uvel')./Uvel');
A4 = Library\Vvel';
Err_Vvel = std(100*(Library*A4 - Vvel')./Uvel');

fprintf('RMS Errors: Speed: %5.2f%%;  Angle: %4.2f Degrees\n', Err_vel, Err_ang);
fprintf('             Uvel: %5.2f%%;    V/U: %4.2f%%\n', Err_Uvel, Err_Vvel);

Vol2S = A1;
Vol2A = A2;
Vol2U = A3;
Vol2V = A4;

if PLOTTING
    %% plot calibration
    FontSize = 14;

    % Find boundaries of data so that we can mask out non-calibrated boundaries
    k = boundary(E1s', E2s');

    pgon1 = polyshape(E1s(k), E2s(k),'Simplify',false);
    pgon  = polybuffer(pgon1, 0.1*(Erange(2)-Erange(1)));  % make a buffer around the voltages

    xlin = linspace(min(pgon.Vertices(:,1)), max(pgon.Vertices(:,1)), 100);
    ylin = linspace(min(pgon.Vertices(:,2)), max(pgon.Vertices(:,2)), 100);
    [X,Y] = meshgrid(xlin, ylin);

    idx = isinterior(pgon, X(:), Y(:));
    idx = reshape(idx, size(X));

    % Plot in terms of Speed and Angle
    figure
    p1 = subplot(2,2,1);
    hold on
    Z = griddata(E1s,E2s,Vel,X,Y,'v4');
    Z(~idx) = nan;
    pcolor(p1, X,Y,Z);
    shading(p1, 'interp')
    colorbar
    plot(E1s, E2s, 'o','MarkerFaceColor','k', 'MarkerEdgeColor', 'k', MarkerSize=5)
    grid on
    title('Speed [m/s]')
    xlabel('E1'); ylabel('E2')
    set(gca,"FontSize", FontSize)
    axis square

    p2 = subplot(2,2,2);
    hold on
    Z = griddata(E1s,E2s,Ang,X,Y,'v4');
    Z(~idx) = nan;
    pcolor(p2, X,Y,Z);
    grid on
    colormap(p2, brewermap([],'PiYG'))
    clim([Angmin Angmax])
    shading('interp')
    colorbar
    plot(E1s, E2s, 'o','MarkerFaceColor','k', 'MarkerEdgeColor', 'k', MarkerSize=5)
    title('Angle [Deg]')
    xlabel('E1'); ylabel('E2')
    set(gca,"FontSize", FontSize)
    axis square

    p3 = subplot(2,2,3);
    hold on
    Z = griddata(E1s,E2s,100*(Library*A1 - Vel')./Vel',X,Y,'v4');
    Z(~idx) = nan;
    pcolor(p3, X,Y,Z)
    grid on
    colormap(p3, brewermap([],'RdBu'))
    shading('interp')
    colorbar
    plot(E1s, E2s, 'o', 'MarkerFaceColor','k', 'MarkerEdgeColor', 'k', MarkerSize=5)
    title('Speed error [rms, %]')
    xlabel('E1'); ylabel('E2')
    set(gca,"FontSize", FontSize)
    axis square

    p4 = subplot(2,2,4);
    hold on
    Z = griddata(E1s,E2s,(Library*A2 - Ang'),X,Y,'v4');
    Z(~idx) = nan;
    pcolor(p4, X,Y,Z)
    grid on
    colormap(p4, brewermap([],'RdBu'))
    shading('interp')
    colorbar
    plot(E1s, E2s, 'o', 'MarkerFaceColor','k', 'MarkerEdgeColor', 'k', MarkerSize=5)
    title('Angle error [rms, deg]')
    xlabel('E1'); ylabel('E2')
    set(gca,"FontSize", FontSize)
    axis square


    % Plot in terms of U and V
    figure
    p1 = subplot(2,2,1);
    hold on
    Z = griddata(E1s,E2s,Uvel,X,Y,'v4');
    Z(~idx) = nan;
    pcolor(p1, X,Y,Z);
    shading(p1, 'interp')
    colorbar
    plot(E1s, E2s, 'o','MarkerFaceColor','k', 'MarkerEdgeColor', 'k', MarkerSize=5)
    grid on
    title('U [m/s]')
    xlabel('E1'); ylabel('E2')
    set(gca,"FontSize", FontSize)
    axis square

    p2 = subplot(2,2,2);
    hold on
    Z = griddata(E1s,E2s,Vvel,X,Y,'v4');
    Z(~idx) = nan;
    pcolor(p2, X,Y,Z);
    grid on
    colormap(p2, brewermap([],'PiYG'))
    shading('interp')
    clim([Vmin Vmax])
    colorbar
    plot(E1s, E2s, 'o','MarkerFaceColor','k', 'MarkerEdgeColor', 'k', MarkerSize=5)
    title('V [m/s]')
    xlabel('E1'); ylabel('E2')
    set(gca,"FontSize", FontSize)
    axis square

    p3 = subplot(2,2,3);
    hold on
    Z = griddata(E1s,E2s,100*(Library*A3 - Uvel')./Uvel',X,Y,'v4');
    Z(~idx) = nan;
    pcolor(p3, X,Y,Z)
    grid on
    colormap(p3, brewermap([],'RdBu'))
    shading('interp')
    colorbar
    plot(E1s, E2s, 'o', 'MarkerFaceColor','k', 'MarkerEdgeColor', 'k', MarkerSize=5)
    title('U-vel error [%]')
    xlabel('E1'); ylabel('E2')
    set(gca,"FontSize", FontSize)
    axis square

    p4 = subplot(2,2,4);
    hold on
    Z = griddata(E1s,E2s,100*(Library*A4 - Vvel')./Uvel',X,Y,'v4');
    Z(~idx) = nan;
    pcolor(p4, X,Y,Z)
    grid on
    colormap(p4, brewermap([],'RdBu'))
    shading('interp')
    colorbar
    plot(E1s, E2s, 'o', 'MarkerFaceColor','k', 'MarkerEdgeColor', 'k', MarkerSize=5)
    title('V/U error [%]')
    xlabel('E1'); ylabel('E2')
    set(gca,"FontSize", FontSize)
    axis square

end



