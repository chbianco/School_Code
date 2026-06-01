% Nathan Wei
% MEAM 6900
% Example code for NREL NWTC-M5 dataset
% Link to data: https://wind.nrel.gov/MetData/135mData/
% Field site information: https://wind.nrel.gov/MetData/Publications/NWTC_135m_MetMasts.pdf
% Metadata information: https://wind.nrel.gov/MetData/Publications/UnofficialGuideToNWTC135mData.pdf 

clear all;

% File location
fileDir = 'C:\Users\YourName\YourDocs\folder';

%% 10-minute averages over the month
load(fullfile(fileDir, '2019_August.mat'));

% Convert UTC datenum to datetime
thyme = all_data.Data_File_Records.date;
time = NaT(size(thyme));
for ti = 1 : length(thyme)
    time(ti) = datetime(floor(thyme(ti)), 'ConvertFrom', 'datenum') ...
        + days(rem(thyme(ti),1)); % convert datenum to datetime object
end
time.TimeZone = 'Etc/UTC'; % add TimeZone field (UTC time)
time.TimeZone = 'America/Denver'; % shift to NREL time zone

% Load variables
U = all_data.Wind_Speed_Horizontal_Sonic_41m.val; % 15m, 41m, 61m, 74m, 100m, 119m

% Plot time series
figure;
plot(time, U);
xlabel('Time');
ylabel('Horizontal Wind Speed (m/s)');

%% 20-Hz data over 10 mins
load(fullfile(fileDir, '08_28_2019_16_00_00_000.mat'));

% Convert UTC datenum to datetime
thyme = time_UTC.val;
time = NaT(size(thyme));
for ti = 1 : length(thyme)
    time(ti) = datetime(floor(thyme(ti)), 'ConvertFrom', 'datenum')...
        + days(rem(thyme(ti),1)); % convert datenum to datetime object
end
time.TimeZone = 'Etc/UTC'; % add TimeZone field (UTC time)
time.TimeZone = 'America/Denver'; % shift to NREL time zone

% Load variables
U = sqrt(Sonic_x_clean_41m.val.^2 + Sonic_y_clean_41m.val.^2); % 15m, 41m, 61m, 74m, 100m, 119m

% Plot time series
figure;
plot(time, U);
xlabel('Time');
ylabel('Horizontal Wind Speed (m/s)');