%% Process x-wire measurements
%
% Kenny Breuer
% April 2025
%
clear all
close all

f_min = 2;  % Minimum frequency for power spectra

% Load raw voltage file which includes the calibrataion data
[Voltage_filename,path]=uigetfile('Xwires_Vol*.mat','select the raw data file');
Voltage_filename=fullfile(path,Voltage_filename);
load(Voltage_filename);

% Number of data points
Npts = length(voltages);

% Compute the calibration data
[Volt2U, Volt2V, Volt2S, Volt2A, Erange] = xwire_compute_calibration(E1_cal, E2_cal, Vel_cal, Ang_cal);

% Size of FFT Spectrum segment
Nseg = 16*round(DigRate/f_min/16);  % Make this a reasonable integer for FFTs

%% Loop through the speeds
figure;

Vels = zeros(Npts, 2, Nspeed_test);
Vel_ave = zeros(2, Nspeed_test);
Vel_std = zeros(2, Nspeed_test);

fprintf('Test data statistics:\n');
for ispeed = 1:Nspeed_test,

    % Convert voltages to velocities:

    E1s = voltages(:,1,ispeed);
    E2s = voltages(:,2,ispeed);

    E1s = E1s';
    E2s = E2s';

    % The calibration library is a fifth order polynomial
    Library = [ ...
        ones(1,Npts); ...
        E1s;    E2s; ...
        E1s.^2; E2s.^2; E1s.*E2s; ...
        E1s.^3; E2s.^3; E1s.*E2s.^2; E1s.^2.*E2s; ...
        E1s.^4; E2s.^4; E1s.*E2s.^3; E1s.^3.*E2s; E1s.^2.*E2s.^2; ...
        E1s.^5; E2s.^5; E1s.*E2s.^4; E1s.^4.*E2s; E1s.^2.*E2s.^3; E1s.^3.*E2s.^2; ...
        ];

    Library = Library';

    % Convert the voltages to velocities using the calibration
    Vels(:,:,ispeed) = [Library*Volt2U Library*Volt2V]/Speed(ispeed);

    % Calculate Mean and Std
    Vel_ave(:,ispeed) = mean(Vels(:,:,ispeed));
    Vel_std(:,ispeed) =  std(Vels(:,:,ispeed));

    % calculate spectrum
    [pxx f] = pwelch(detrend(Vels(:,:,ispeed)), Nseg, Nseg/2, Nseg, DigRate);

    % report the results
    fprintf('Speed: %d/%d: Uo: %6.2f; U/Uo: %5.2f [+/- %5.3f]; V/Uo: %5.3f [+/- %5.3f]; \n', ...
        ispeed, Nspeed_test, Speed(ispeed), Vel_ave(1,ispeed), Vel_std(1,ispeed), Vel_ave(2,ispeed), Vel_std(2,ispeed));

    % Plot results
    iplot =(ispeed-1)*3;
    FontSize = 14;

    % Plot the spectra
    subplot(Nspeed_test,3,iplot+1)
    semilogx(f,20*log10(pxx))
    hold on
    xlabel('Frequency, [Hz]')
    ylabel('Power [dB]')
    title(sprintf('%6.2f m/s', Speed(ispeed)));
    grid on
    legend('u''', 'v''', 'Location', 'SouthWest')
    set(gca,'FontSize', FontSize)

    % Plot the distribution of velocity fluctuations along with the calibration data
    subplot(Nspeed_test,3,iplot+2)
    plot((decimate(Vels(:,1,ispeed),10)-Vel_ave(1,ispeed)), ...
         (decimate(Vels(:,2,ispeed),10)-Vel_ave(2,ispeed)), 'bo',MarkerSize=1);
    xlabel('u''/U_o');
    ylabel('v''/U_o');
    axis equal
    xline(0)
    yline(0)
    set(gca,'FontSize', FontSize)

    % Plot the distribution of voltages, compared with the calibration
    % data  (dont plot all the data)
    subplot(Nspeed_test,3,iplot+3)
    plot(decimate(voltages(:,1,ispeed),10), decimate(voltages(:,2,ispeed),10), 'bo',MarkerSize=1);
    hold on
    plot(reshape(E1_cal,1,Nang_cal*Nspeed_cal), reshape(E2_cal,1,Nang_cal*Nspeed_cal),'o', ...
        'MarkerFaceColor','r', 'MarkerEdgeColor', 'r', MarkerSize=5);
    xlabel('E1');
    ylabel('E2');
    axis square
    set(gca,'FontSize', FontSize)
end

% Save velocities and experimental data
time_string = datetime('now','TimeZone','local','Format','d-MMM-y-HH-mm-ss');
velocity_filename = sprintf('Xwires_Vel_%s',time_string);
if ~exist('calib_filename', 'var')
    calib_filename = 'unknown';
end
eval(['save ' velocity_filename ' Nspeed_test voltages Vels Speed Temp DigRate Density Pressure Viscosity' ...
            ' E1_cal E2_cal Vel_cal Ang_cal Temp_cal Nspeed_cal Nang_cal calib_filename' ]);
fprintf('File saved: %s\n', velocity_filename);
