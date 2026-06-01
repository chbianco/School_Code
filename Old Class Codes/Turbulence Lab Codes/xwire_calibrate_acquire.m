%% Code to measure hot wire voltages at different speeds and angles
%
% Kenny Breuer, Jan 2025
%
% Make sure the Thor stage is clean (wipe with alcohol q-tip)
%
% NOTE that the hot wires have high voltage at low speed and vice versa,
% and so the code flips that, just to make it more intuitive.  For the DAQ,
% the voltages need to be between +/- 5V
%
% At low speed the zero angle voltage should be about 3-4V

% 
% Pressures in the ELD tunnel should be as low as possible to get a low speed, 
% up to about 40 PSI, for a top calibration speed of ~30 m/s

clear all

%% Parameters for the acquisition
Total_Time = 5;   % Total time at each speed to take data
DigRate = 10240;   % Digitizing frquency (Hz)

Nspeed_cal = 5;
Nangles = 7;
Angle_Max = 30;

NI_DEV = 'cDAQ1Mod1';

AFAM_TUNNEL = 0;

fprintf('Each speed will take %3.1f minutes\n', Total_Time*Nangles/60);

if AFAM_TUNNEL
    input_string = 'RPM';
else
    input_string = 'pressure [PSI]';
    Agilent_State;  % Initialize Agilent
end

disp(['Enter ' input_string ' for wind tunnel speeds'])
Set_min = input(['Enter first ' input_string ': ']);
Set_max = input(['Enter last ' input_string ': ']);
RPM = round(linspace(Set_min, Set_max, Nspeed_cal));

% Speeds must be >= 6 for a fifth order polynomial calibration
Angles = linspace(-Angle_Max, Angle_Max, Nangles);

%% Set up the NI-DAQ
s = daq('ni');
s.Rate = DigRate;
addinput(s, NI_DEV, 'ai0', 'Voltage');
addinput(s, NI_DEV, 'ai1', 'Voltage');

% Zero out the arrays
E_ave    = zeros(2,Nangles, Nspeed_cal);
E_std    = zeros(2,Nangles, Nspeed_cal);
Vel      = zeros(  Nangles, Nspeed_cal);
Temp     = zeros(  Nangles, Nspeed_cal);

% Loop through the speeds

disp('Is the probe at zero degrees to the flow? - hit any key to confirm')
pause;
THOR = Thor_Rotation_Stage('init');  % Initialize stage
Theta0 = THOR.posDegree;  % get zero position

figure
axis([-6 6 -6 6])
axis equal

for ispeed = 1:Nspeed_cal,

    if AFAM_TUNNEL
        fprintf('Setting RPM to %d...', Set(ispeed));
        VFD_set_RPM(Set(ispeed))
        pause(20);
    else
        fprintf('Speed %d/%d; set pressure to %d.  [Hit any key when ready to continue]', ...
            ispeed, Nspeed_cal, RPM(ispeed))
        pause;
        fprintf('\n')
    end

    for iangle = 1:Nangles
        THOR = Thor_Rotation_Stage('gp'); % get current position
        Current_Theta = THOR.posDegree - Theta0;
        angle_to_move  = Angles(iangle) - Current_Theta;

        THOR = Thor_Rotation_Stage('mr', angle_to_move);
        pause(1);

        % Measure ADC
        data = read(s, seconds(Total_Time));

        % get true position
        THOR = Thor_Rotation_Stage('gp'); % get current position
        Theta(iangle,ispeed) = THOR.posDegree - Theta0;


        % Get the speed and temperature from the AFAM software
        if AFAM_TUNNEL
            Vel(iangle,ispeed)  = AFAM_Tunnel.Speed;      % Speed from Pitot Tube
            Temp(iangle,ispeed) = AFAM_Tunnel.TestTemp;   % Temperature of test section
        else
            [Temp(iangle,ispeed) P_diff P_abs Density Viscosity Vel(iangle,ispeed)] = Agilent_State;
        end

        % Calculate Mean and Std
        % Change the sign so that the signal increases with increasing
        % velocity
        E_ave(:,iangle, ispeed) = -mean(table2array(data));
        E_std(:,iangle, ispeed) =  std(table2array(data));

        % report the results
        fprintf(['   V: %' ...
            '5.2f; T: %5.2f; Ang: %6.2f: Ch 1: %5.2f [+/- %6.4f]; Ch2: %5.2f [+/- %6.4f]\n'], ...
            Vel(iangle,ispeed), Temp(iangle, ispeed), Theta(iangle,ispeed), E_ave(1, iangle, ispeed), E_std(1, iangle, ispeed), ...
            E_ave(2, iangle, ispeed), E_std(2, iangle, ispeed));

        plot(E_ave(1, iangle, ispeed), E_ave(2, iangle, ispeed), 'bo', 'MarkerSize', 4, 'MarkerFaceColor','b')
        xlabel('Channel -E1')
        ylabel('Channel -E2')
        hold on
        set(gca, 'FontSize', 16)
        axis 'square'

    end
end

% Move back to zero
THOR = Thor_Rotation_Stage('gp'); % get current position
Current_Theta = THOR.posDegree - Theta0;
Thor_Rotation_Stage('mr',-Current_Theta);

if AFAM_TUNNEL
    VFD_stop;
end

Thor_Rotation_Stage('close');

%% Save the variables for later processing
E1_cal = E_ave(1,:,:);
E2_cal = E_ave(2,:,:);
Ang_cal = Theta;
Vel_cal = Vel;
Temp_cal = Temp;
Nang_cal = Nangles;

%% Save data
time_string = datetime('now','TimeZone','local','Format','d-MMM-y-HH-mm-ss');
calib_filename = sprintf('Xwires_Cal_%s',time_string);
eval(['save ' calib_filename ' E1_cal E2_cal Vel_cal Ang_cal Temp_cal Nspeed_cal Nang_cal calib_filename']);
fprintf('File saved: %s\n', calib_filename);






