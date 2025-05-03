%% Measure x-wire wire voltages in wind tunnel
%
% Kenny Breuer
% Jan 2024
%
% expanded in April 2025 for X-wires
%
clear all
close all


%% Parameters
DigRate = 51200; % Hz - the dDAQ has fixed digitizing rate
Npts = 500000;
Nspeed_test= 3;
f_min = 2; % Minimum frequency to resolve

Nch = 2;
NI_DEV = 'cDAQ1Mod1';
Channel0 = 'ai0';
Channel1 = 'ai1';


% Load the raw calibration file
uiopen('Xwires_Cal*.mat')

%% Get on with the experiment

E1 = E1_cal(1,:,:);
E2 = E2_cal(1,:,:);

%reduce the rank for data fitting
E1s = reshape(E1(1,:,:),1, Nang_cal*Nspeed_cal);
E2s = reshape(E2(1,:,:),1, Nang_cal*Nspeed_cal);
Vel = reshape(Vel_cal,  1, Nang_cal*Nspeed_cal);
Ang = reshape(Ang_cal,  1, Nang_cal*Nspeed_cal);

%fprintf('Calibration data from %s\n', [pathname filename])
fprintf('Velocity: %5.1f - %5.1f; Range: %5.2f [m/s]\n', min(Vel), max(Vel), max(Vel) - min(Vel));
fprintf('   Angle: %5.1f - %5.1f; Range: %5.2f [deg]\n', min(Ang), max(Ang), max(Ang) - min(Ang));
fprintf('      E1: %5.2f - %5.2f; Range: %5.2f [V]\n',   min(E1s), max(E1s), max(E1s) - min(E1s));
fprintf('      E2: %5.2f - %5.2f; Range: %5.2f [V]\n',   min(E2s), max(E2s), max(E2s) - min(E2s));

%% Set up the NI-DAQ
s = daq.createSession('ni');
s.NumberOfScans = Npts;
s.Rate = DigRate;
addAnalogInputChannel(s, NI_DEV, Channel0, 'Voltage');  
addAnalogInputChannel(s, NI_DEV, Channel1, 'Voltage'); 

% Size of FFT Spectrum segment
Nseg = 16*round(DigRate/f_min/16);  % Make this a reasonable integer for FFTs

%% Initialize the Agilent box
P_zero = 45; %input('Set Velocity to zero; enter tare voltage on Channel 201 (Diff Press) in mV: ');    
Agilent_State(P_zero);

fprintf('Each speed will take approx %d seconds\n', round(Npts/DigRate));

% Zero out the arrays
clear V_ave V_std vel legend_text;
voltages = zeros(Npts,2,Nspeed_test);

%% Loop through the speeds
figure;

for ispeed = 1:Nspeed_test,
    fprintf('Speed: %d/%d:  Set tunnel speed; wait for speed to settle, then hit any key to start ADC', ...
        ispeed, Nspeed_test);
    pause
    fprintf(' .. digitizing ..\n')
    pause(0.5)
    data = startForeground(s);

    % Flip the voltages to be consistent with the calibration
    data = -data;

    [Temp(ispeed), P_diff(ispeed), Pressure(ispeed), Density(ispeed), Viscosity(ispeed), Speed(ispeed)] = ...
        Agilent_State;
    
    fprintf('   Vel: %6.2f [m/s]; Temp: %5.2f [C]\n',Speed(ispeed), Temp(ispeed));
    for ich = 1:Nch,
        % Save the data for later
        voltages(:,ich, ispeed) = data(:,ich);

        % Calculate Mean and Std
        V_ave(ich,ispeed) = mean(data(:,ich));
        V_std(ich,ispeed) =  std(data(:,ich));

        % calculate spectrum
        [pxx f] = pwelch(detrend(data(:,ich)), Nseg, Nseg/2, Nseg, DigRate);

        % report the results
        fprintf('   Channel: %d; Mean-V: %6.3f [+/- %7.4f]; Range: %5.2f - %5.2f [V].\n', ...
            ich, V_ave(ich,ispeed), V_std(ich,ispeed), min(data(:,ich)), max(data(:,ich)));

        % Plot results
        iplot =(ispeed-1)*2;

        subplot(Nspeed_test,2,iplot+1)
        semilogx(f,10*log10(pxx))
        hold on
        xlabel('Frequency, [Hz]')
        ylabel('Voltage Power [dB]')
        title(sprintf('%6.2f m/s', Speed(ispeed)));

        % Plot the distribution of voltages along with the calibration data
        subplot(Nspeed_test,2,iplot+2)
        plot(data(:,1), data(:,2), 'ko');
        hold on
        plot(E1s, E2s, 'o', 'MarkerFaceColor','r', 'MarkerEdgeColor', 'r', MarkerSize=2);
        xlabel('E1');
        ylabel('E2');
    end
end

% Clear the data acq data structure - doesnt save well
clear s

%% Save data including the calibration data
time_string = datetime('now','TimeZone','local','Format','d-MMM-y-HH-mm-ss');
filename = sprintf('Xwires_Vol_%s',time_string);
% eval(['save ' filename]);
eval(['save ' filename ' Nspeed_test voltages Speed Temp DigRate Density Pressure Viscosity' ...
            ' E1_cal E2_cal Vel_cal Ang_cal Temp_cal Nspeed_cal Nang_cal calib_filename' ]);
fprintf('File saved: %s\n', filename);
