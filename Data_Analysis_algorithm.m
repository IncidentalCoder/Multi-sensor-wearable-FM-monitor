% Code for analyzing fetal movement data obtained from the new multi-sensor
% fetal movement monitor
%

%% DATA LOADING ===========================================================
% This part of the code loads the data files for processing. Loading of
% single or multiple data files are allowed. When a data file is loaded, a
% variable named data_var is loaded in the workspace from where the sensor
% and sensation data for that particular session is extracted and stored in
% a cell matrix.

clear; clc;

% Starting notification
disp('Data file loading is going on...');

% Adding path for function files
addpath(genpath('Function_files'))

% Locating the data files
[data_file_names, path_name] = uigetfile('*.mat', 'Select the data files','MultiSelect', 'on'); % Returns file names with the extension in a cell, and path name
addpath(path_name); % Adds the path to the file

% Counting the number of data files to be loaded
%   data_file_names is a string in case of a single data file selected.
%   Otherwise, it is a cell variable containing 1 namein each cell    
if ischar(data_file_names) 
    n_data_files = 1;    
else
    n_data_files = length(data_file_names); 
end

% Variable decleration
sensation_data_SD1 = cell(1,n_data_files); 
sensation_data_SD2 = cell(1,n_data_files); 

Acstc_data1 = cell(1,n_data_files);  Acstc_data2 = cell(1,n_data_files);
Aclm_data1 = cell(1,n_data_files);   Aclm_data2 = cell(1,n_data_files);
Pzplt_data1 = cell(1,n_data_files);  Pzplt_data2 = cell(1,n_data_files);
Flexi_data = cell(1,n_data_files);   IMU_data = cell(1,n_data_files);

% Known parameters
Fs_sensor = 1024; % Frequency of sensor data sampling in Hz
Fs_sensation = 1024; % Frequency of sensation data sampling in Hz
n_sensor = 8; % Total number of sensors

% Loading the data files
for i = 1 : n_data_files
    
    if n_data_files == 1
        load(data_file_names); % A cell matrix named data_var will be loaded, which contains all the data 
    else
        load(data_file_names{i}); % A cell matrix named data_var will be loaded, which contains all the data
    end
    
    sensor_data_SD1 = data_var{1,1}; 
    sensor_data_SD2 = data_var{1,3};
    sensation_data_SD1{i} = data_var{1,2};
    sensation_data_SD2{i} = data_var{1,4};
    
    % Extracting sensor data
    Flexi_data{i} = sensor_data_SD1(:, 1);  % Force sensor data to measure tightness of the belt
    Pzplt_data1{i} = sensor_data_SD1(:, 2); % Left Piezo-plate data to measure abdomenal vibration
    Pzplt_data2{i} = sensor_data_SD1(:, 3); % Right Piezo-plate data to measure abdomenal vibration
    Acstc_data1{i} = sensor_data_SD1(:, 4); % Left Acoustic sensor data to measure abdomenal vibration
    IMU_data{i} = sensor_data_SD1(:, 5:7);  % Accelerometer data to measure maternal body movement
    
    Acstc_data2{i} = sensor_data_SD2(:, 1); % Right Acoustic sensor data to measure abdomenal vibration
    Aclm_data1{i} = sensor_data_SD2(:, 5:7);% Left Accelerometer data to measure abdomenal vibration
    Aclm_data2{i} = sensor_data_SD2(:, 2:4);% Right Accelerometer data to measure abdomenal vibration
    
    % Resultant of acceleration
    IMU_data{i} = sqrt(sum(IMU_data{i}.^2,2));
    Aclm_data1{i} = sqrt(sum(Aclm_data1{i}.^2,2));
    Aclm_data2{i} = sqrt(sum(Aclm_data2{i}.^2,2));
    
end

fprintf('In total, %.0f data files have been uploaded.\n\n', n_data_files);
%

%% ============================= PRE-PROCESSING ===========================
% Pre-processing steps includes filtering/detrending and trimming the data 
% Triming is done to remove some initial and final unwanted data
% Force data analysis is also included in this step

% --------------------------- Variable decleration -----------------------%
Acstc_data1_fltd = cell(1,n_data_files); Acstc_data2_fltd = cell(1,n_data_files);
Aclm_data1_fltd = cell(1,n_data_files);  Aclm_data2_fltd = cell(1,n_data_files);
Pzplt_data1_fltd = cell(1,n_data_files); Pzplt_data2_fltd = cell(1,n_data_files);
Flexi_data_fltd = cell(1,n_data_files);  IMU_data_fltd = cell(1,n_data_files);

sensation_data_SD1_trimd = cell(1,n_data_files); 
sensation_data_SD2_trimd = cell(1,n_data_files); 

Acstc_data1_trimd = cell(1,n_data_files); Acstc_data2_trimd = cell(1,n_data_files);
Aclm_data1_trimd  = cell(1,n_data_files); Aclm_data2_trimd  = cell(1,n_data_files);
Pzplt_data1_trimd = cell(1,n_data_files); Pzplt_data2_trimd = cell(1,n_data_files);
Flexi_data_trimd  = cell(1,n_data_files); IMU_data_trimd    = cell(1,n_data_files);

% ---------------------------- Filter design -----------------------------%
%   3 types of filters are designed- bandpass filter, low-pass filter, and IIR notch filter.

filter_order = 10;

%   Bandpass filter
%       A band-pass filter with a passband of 1-20 Hz is disigned for the fetal fetal movement data
%       Another band-pass filer with a passband of 1-10 Hz is designed for the IMU data
lowCutoff_FM = 1;
highCutoff_FM = 30;
lowCutoff_IMU = 1;
highCutoff_IMU = 10;

%       Transfer function-based desing
[b_FM,a_FM] = butter(filter_order/2,[lowCutoff_FM highCutoff_FM]/(Fs_sensor/2),'bandpass');
[b_IMU,a_IMU] = butter(filter_order/2,[lowCutoff_IMU highCutoff_IMU]/(Fs_sensor/2),'bandpass');

%       Zero-Pole-Gain-based design
[z_FM,p_FM,k_FM] = butter(filter_order/2,[lowCutoff_FM highCutoff_FM]/(Fs_sensor/2),'bandpass'); % filter order for bandpass filter is twice the value of 1st parameter
[sos_FM,g_FM] = zp2sos(z_FM,p_FM,k_FM); % Convert zero-pole-gain filter parameters to second-order sections form
 
[z_IMU,p_IMU,k_IMU] = butter(filter_order/2,[lowCutoff_IMU highCutoff_IMU]/(Fs_sensor/2),'bandpass');
[sos_IMU,g_IMU] = zp2sos(z_IMU,p_IMU,k_IMU);
%

% 	Low-pass filter 
%       This filter is used for the force sensor data only
highCutoff_force = 10;
[z_force,p_force,k_force] = butter(filter_order,highCutoff_force/(Fs_sensor/2),'low'); % 
[sos_force,g_force] = zp2sos(z_force,p_force,k_force); % Convert zero-pole-gain filter parameters to second-order sections form

% ------------- Loop for filtering and trimming the data -----------------%
% Trim settings
start_removal_period = 30; % Removal period in second
end_removal_period = 30; % Removal period in second

% Starting notification
disp('Pre-processing is going on with the following setting-');
fprintf(['\tFilter order: %.1f \n\tIMU band-pass: %.1f-%.1f Hz \n\tFM band-pass: %.1f-%.1f Hz\n\tForce sensor low-pass: %.1f Hz\n ' ...
    'Data trimming at the start: %.1f \n Data trimming at the end: %.1f ...\n'], filter_order, lowCutoff_IMU, highCutoff_IMU, ...
    lowCutoff_FM, highCutoff_FM, highCutoff_force, start_removal_period, end_removal_period );

for i = 1 : n_data_files
    
    % ------------------------- Data filtering ---------------------------%
    Acstc_data1_fltd{i} = Acstc_data1{i};
    Acstc_data2_fltd{i} = Acstc_data2{i}; 
    Aclm_data1_fltd{i}  = Aclm_data1{i};
    Aclm_data2_fltd{i}  = Aclm_data2{i};
    Pzplt_data1_fltd{i} = Pzplt_data1{i};
    Pzplt_data2_fltd{i} = Pzplt_data2{i};
    Flexi_data_fltd{i}  = Flexi_data{i};
    IMU_data_fltd{i}    = IMU_data{i};
    %

    % Low-pass filtering
    Flexi_data_fltd{i}  = filtfilt(sos_force, g_force, Flexi_data_fltd{i});
    
    % Bandpass filtering
    Acstc_data1_fltd{i} = filtfilt(sos_FM, g_FM, Acstc_data1_fltd{i}); % Zero-phase filtering, and filter oder is twice than the designed sos
    Acstc_data2_fltd{i} = filtfilt(sos_FM, g_FM, Acstc_data2_fltd{i}); 
    Aclm_data1_fltd{i}  = filtfilt(sos_FM, g_FM, Aclm_data1_fltd{i});
    Aclm_data2_fltd{i}  = filtfilt(sos_FM, g_FM, Aclm_data2_fltd{i});
    Pzplt_data1_fltd{i} = filtfilt(sos_FM, g_FM, Pzplt_data1_fltd{i});
    Pzplt_data2_fltd{i} = filtfilt(sos_FM, g_FM, Pzplt_data2_fltd{i});    
    IMU_data_fltd{i}    = filtfilt(sos_IMU,g_IMU,IMU_data_fltd{i});
    %
   
    % ---------------------------- Data trimming -------------------------%    
    % Trimming of raw data    
    sensation_data_SD1_trimd{i} = sensation_data_SD1{i}((start_removal_period*Fs_sensation + 1) : (end-end_removal_period*Fs_sensation));
    sensation_data_SD2_trimd{i} = sensation_data_SD2{i}((start_removal_period*Fs_sensation + 1) : (end-end_removal_period*Fs_sensation));
    
    Pzplt_data1_trimd{i} = Pzplt_data1{i}((start_removal_period*Fs_sensor + 1) : (end-end_removal_period*Fs_sensor));
    Pzplt_data2_trimd{i} = Pzplt_data2{i}((start_removal_period*Fs_sensor + 1) : (end-end_removal_period*Fs_sensor));
    Acstc_data1_trimd{i} = Acstc_data1{i}((start_removal_period*Fs_sensor + 1) : (end-end_removal_period*Fs_sensor));
    Acstc_data2_trimd{i} = Acstc_data2{i}((start_removal_period*Fs_sensor + 1) : (end-end_removal_period*Fs_sensor));
    Aclm_data1_trimd{i}  = Aclm_data1{i}((start_removal_period*Fs_sensor + 1)  : (end-end_removal_period*Fs_sensor));
    Aclm_data2_trimd{i}  = Aclm_data2{i}((start_removal_period*Fs_sensor + 1)  : (end-end_removal_period*Fs_sensor));
    Flexi_data_trimd{i}  = Flexi_data{i}((start_removal_period*Fs_sensor + 1)  : (end-end_removal_period*Fs_sensor));
    IMU_data_trimd{i}    = IMU_data{i}((start_removal_period*Fs_sensor + 1)    : (end-end_removal_period*Fs_sensor));
    %
    
    % Trimming of filtered data
    Pzplt_data1_fltd{i} = Pzplt_data1_fltd{i}((start_removal_period*Fs_sensor + 1) : (end-end_removal_period*Fs_sensor));
    Pzplt_data2_fltd{i} = Pzplt_data2_fltd{i}((start_removal_period*Fs_sensor + 1) : (end-end_removal_period*Fs_sensor));
    Acstc_data1_fltd{i} = Acstc_data1_fltd{i}((start_removal_period*Fs_sensor + 1) : (end-end_removal_period*Fs_sensor));
    Acstc_data2_fltd{i} = Acstc_data2_fltd{i}((start_removal_period*Fs_sensor + 1) : (end-end_removal_period*Fs_sensor));
    Aclm_data1_fltd{i}  = Aclm_data1_fltd{i}((start_removal_period*Fs_sensor + 1)  : (end-end_removal_period*Fs_sensor));
    Aclm_data2_fltd{i}  = Aclm_data2_fltd{i}((start_removal_period*Fs_sensor + 1)  : (end-end_removal_period*Fs_sensor));
    Flexi_data_fltd{i}  = Flexi_data_fltd{i}((start_removal_period*Fs_sensor + 1)  : (end-end_removal_period*Fs_sensor));
    IMU_data_fltd{i}    = IMU_data_fltd{i}((start_removal_period*Fs_sensor + 1)    : (end-end_removal_period*Fs_sensor));
    %
    
    % Equalizing the length of SD1 and SD2 data sets
    min_length_sensation = min(length(sensation_data_SD1_trimd{i}),length(sensation_data_SD2_trimd{i}));
    min_length_sensor = min(length(Acstc_data1_fltd{i}),length(Acstc_data2_fltd{i}));
    % min_length_sensation and min_length_sensor are equal except for the case of data sets taken with DAQ version 1.0 
    
    sensation_data_SD1_trimd{i} = sensation_data_SD1_trimd{i}(1:min_length_sensation);
    sensation_data_SD2_trimd{i} = sensation_data_SD2_trimd{i}(1:min_length_sensation);
    
    Pzplt_data1_trimd{i} = Pzplt_data1_trimd{i}(1:min_length_sensor);
    Pzplt_data2_trimd{i} = Pzplt_data2_trimd{i}(1:min_length_sensor);
    Acstc_data1_trimd{i} = Acstc_data1_trimd{i}(1:min_length_sensor);
    Acstc_data2_trimd{i} = Acstc_data2_trimd{i}(1:min_length_sensor);
    Aclm_data1_trimd{i}  = Aclm_data1_trimd{i}(1:min_length_sensor);
    Aclm_data2_trimd{i}  = Aclm_data2_trimd{i}(1:min_length_sensor);
    Flexi_data_trimd{i}  = Flexi_data_trimd{i}(1:min_length_sensor);
    IMU_data_trimd{i}    = IMU_data_trimd{i}(1:min_length_sensor);
    
    Pzplt_data1_fltd{i} = Pzplt_data1_fltd{i}(1:min_length_sensor);
    Pzplt_data2_fltd{i} = Pzplt_data2_fltd{i}(1:min_length_sensor);
    Acstc_data1_fltd{i} = Acstc_data1_fltd{i}(1:min_length_sensor);
    Acstc_data2_fltd{i} = Acstc_data2_fltd{i}(1:min_length_sensor);
    Aclm_data1_fltd{i}  = Aclm_data1_fltd{i}(1:min_length_sensor);
    Aclm_data2_fltd{i}  = Aclm_data2_fltd{i}(1:min_length_sensor);
    Flexi_data_fltd{i}  = Flexi_data_fltd{i}(1:min_length_sensor);
    IMU_data_fltd{i}    = IMU_data_fltd{i}(1:min_length_sensor);
    %     
end
%

% ----------------------- Extracting general info -------------------------
% Extracting force sensor data during each recording session
Force_data_sample = cell(1,n_data_files);
Force_mean = zeros(n_data_files, 1); 
Force_signal_power = zeros(1, n_data_files);
% sample_size = 30; % sample size in seconds

for i = 1 : n_data_files    
    %Force_data_sample{i} = Flexi_data_fltd{i}(1:sample_size*Fs_sensor);
    Force_data_sample{i} = abs(Flexi_data_fltd{i});
    Force_mean(i) = mean(Force_data_sample{i});
    Force_signal_power(i) = sum(Force_data_sample{i}.^2)/length(Force_data_sample{i});
end
%

% Duration of each data sample
duration_data_files = zeros(n_data_files,1);
for i = 1:n_data_files
    duration_data_files(i)= length(Acstc_data1{i})/(Fs_sensor*60); % Duration in minutes
end
%

% Ending notification 
fprintf('Data filtering and trimming are completed.\n\n')
%

%% PERFORMANCE EVALUATION BASED ON OPTIMUM MULTIPLIER =====================
% Starting notification
disp('Performance analysis is going on... ')

% Parameters for FM data segmentation
n_FM_sensors = 6; % number of FM sensors
FM_dilation_time = 3; % Dilation time for FM signal in second
FM_min_SN = [30, 30, 60, 60, 50, 50];

% Parameters for creating sensation map and detection matching
ext_backward = 5.0; % Backward extension length in second
ext_forward = 2.0 ; % Forward extension length in second

% Variable decleration
IMU_map = cell(1,n_data_files);   
M_sntn_map = cell(1,n_data_files);  
n_activity = zeros(n_data_files, 1); % Cell matrix to hold the number of fetal activities
n_Maternal_detected_movement = zeros(n_data_files, 1);
threshold = zeros(n_data_files, n_FM_sensors); % Variable for thresold

TPD_all_indv_sensors_indv = cell(1, n_FM_sensors);
FPD_all_indv_sensors_indv = cell(1, n_FM_sensors);
TND_all_indv_sensors_indv = cell(1, n_FM_sensors);
FND_all_indv_sensors_indv = cell(1, n_FM_sensors);

for i = 1 : n_data_files

    % Starting notification
    fprintf('Current data file: %.0f/%.0f\n', i, n_data_files)

    % --------- Segmentaiton of IMU data and creation of IMU_map ---------%
    % get_IMU_map() function is used here. It segments and dilates the IMU 
    % data and returns the resultant data as IMU map. Settings for the
    % segmentation and dilation are given inside the function.
    %   Input variables:  IMU_data- A column vector of IMU data 
    %                     data_file_names- a char variable with data file name
    %   Output variables: IMU_map- A column vector with segmented IMU data 

    if n_data_files == 1
        IMU_map{i} = get_IMU_map(IMU_data_fltd{i}, data_file_names, Fs_sensor);
    else
        IMU_map{i} = get_IMU_map(IMU_data_fltd{i}, data_file_names{i}, Fs_sensor);
    end

    % ----------------------- Creation of M_sensation_map ----------------%
    % get_sensation_map() function is used here. This function gives the
    % sensation map by dilating every detection to past and future. It also 
    % revomes the windows that overlaps with IMU_map. Settings for the 
    % extension are given inside the function.
    %   Input variables-  sensation_data- a column vector with the sensation data
    %                     IMU_map- a column vector
    %                     ext_backward, ext_forward- scalar values
    %                     Fs_sensor,Fs_sensation- scalar values    %                     
    %   Output variables: M_sntn_map- a column vector with all the sensation windows joined together 

    M_sntn_map{i} = get_sensation_map(sensation_data_SD1_trimd{i}, IMU_map{i}, ext_backward, ext_forward, Fs_sensor, Fs_sensation);

    % ---------------- Determination of detection statistics --------------
    M_sntn_Map_labeled = bwlabel(M_sntn_map{i});
    n_activity(i) = length(unique(M_sntn_Map_labeled)) - 1; % 1 is deducted to remove the first element, which is 0

    sensation_data_labeled = bwlabel(sensation_data_SD1_trimd{i});
    n_Maternal_detected_movement(i) = length(unique(sensation_data_labeled)) - 1; % 1 is deducted to remove the initial value of 0
    %

    % ~~~~~~~~~~~~~~~~~~~~~~ Segmentation of FM data ~~~~~~~~~~~~~~~~~~~~~%
    % get_segmented_data() function will be used here, which will
    % threshold the data, remove body movement, and dilate the data.
    % Setting for the threshold and dilation are given in the function.
    %   Input variables:  sensor_data- a cell variable;
    %                     min_SN- a vector/ a scalar;
    %                     IMU_map- a column vector
    %                     Fs_sensor, FM_dilation_time- a scalar
    %   Output variables: sensor_data_sgmntd- a cell variable of same size as the input variable sensor_data_fltd;
    %                     h- a vector

    sensor_data_fltd = {Aclm_data1_fltd{i}, Aclm_data2_fltd{i}, Acstc_data1_fltd{i}, Acstc_data2_fltd{i},  ...
        Pzplt_data1_fltd{i}, Pzplt_data2_fltd{i}};
    [sensor_data_sgmntd, threshold(i,:)] = get_segmented_data(sensor_data_fltd, FM_min_SN, IMU_map{i}, FM_dilation_time, Fs_sensor);

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~ Sensor fusion ~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    % First fuse the left and right sensors of each type using OR
    sensor_data_sgmntd_Aclm = {double(sensor_data_sgmntd{1} | sensor_data_sgmntd{2})};
    sensor_data_sgmntd_Acstc = {double(sensor_data_sgmntd{3} | sensor_data_sgmntd{4})};    
    sensor_data_sgmntd_Pzplt = {double(sensor_data_sgmntd{5} | sensor_data_sgmntd{6})};

    % Then apply different sensor fusion scheme to combine different types of sensors
    % SENSOR FUSION SCHEME 1: Combination based on logical OR operation
    %   In this scheme, data fusion is performed before dilation.
    %   Combined data are stored as a cell to make it compatable with the
    %   function related to matcing with maternal sensation
    sensor_data_sgmntd_Aclm_OR_Acstc  = {double(sensor_data_sgmntd_Aclm{1} | sensor_data_sgmntd_Acstc{1})};
    sensor_data_sgmntd_Aclm_OR_Pzplt  = {double(sensor_data_sgmntd_Aclm{1} | sensor_data_sgmntd_Pzplt{1})};
    sensor_data_sgmntd_Acstc_OR_Pzplt = {double(sensor_data_sgmntd_Acstc{1} | sensor_data_sgmntd_Pzplt{1})};   

    sensor_data_sgmntd_cmbd_all_sensors_OR = {double(sensor_data_sgmntd_Aclm{1} | sensor_data_sgmntd_Acstc{1} | sensor_data_sgmntd_Pzplt{1})};
    %

    % SENSOR FUSION SCHEME 2: Combination based on logical AND operation
    %   In this scheme, data fusion is performed after dilation.
    %   Combined data are stored as a cell to make it compatable with the
    %   function related to matcing with maternal sensation

    sensor_data_sgmntd_Aclm_AND_Acstc  = {double(sensor_data_sgmntd_Aclm{1} & sensor_data_sgmntd_Acstc{1})};
    sensor_data_sgmntd_Aclm_AND_Pzplt  = {double(sensor_data_sgmntd_Aclm{1} & sensor_data_sgmntd_Pzplt{1})};
    sensor_data_sgmntd_Acstc_AND_Pzplt = {double(sensor_data_sgmntd_Acstc{1} & sensor_data_sgmntd_Pzplt{1})}; 

    sensor_data_sgmntd_cmbd_all_sensors_AND = {double(sensor_data_sgmntd_Aclm{1} & sensor_data_sgmntd_Acstc{1} & sensor_data_sgmntd_Pzplt{1})};
    %

    % SENSOR FUSION SCHEME 3: Combination based on detection by atleaset n sensors
    %   In this scheme, data fusion is performed after dilation.
    %   All the sensor data are first combined with logical OR. Each
    %   non-zero sengment in the combined data is then checked against
    %   individual sensor data to find its presence in that data set.
    %   Combined data are stored as a cell to make it compatable with the
    %   function related to matcing with maternal sensation.
    sensor_data_sgmntd_cmbd_all_sensors = {double(sensor_data_sgmntd{1} | sensor_data_sgmntd{2} | sensor_data_sgmntd{3} | ...
        sensor_data_sgmntd{4} | sensor_data_sgmntd{5} | sensor_data_sgmntd{6})};
    sensor_data_sgmntd_cmbd_all_sensors_labeled = bwlabel(sensor_data_sgmntd_cmbd_all_sensors{1});
    n_label = length(unique(sensor_data_sgmntd_cmbd_all_sensors_labeled)) - 1; % Number of labels in the sensor_data_cmbd_all_sensors_labeled

    sensor_data_sgmntd_cmbd_multi_type_sensors_OR = [sensor_data_sgmntd_Acstc, sensor_data_sgmntd_Aclm, sensor_data_sgmntd_Pzplt];

    % Initialization of variables
    sensor_data_sgmntd_atleast_1_sensor{1} = zeros(length(sensor_data_sgmntd{1}),1); 
    sensor_data_sgmntd_atleast_2_sensor{1} = zeros(length(sensor_data_sgmntd{1}),1);
    sensor_data_sgmntd_atleast_3_sensor{1} = zeros(length(sensor_data_sgmntd{1}),1);
    sensor_data_sgmntd_atleast_4_sensor{1} = zeros(length(sensor_data_sgmntd{1}),1);
    sensor_data_sgmntd_atleast_5_sensor{1} = zeros(length(sensor_data_sgmntd{1}),1);
    sensor_data_sgmntd_atleast_6_sensor{1} = zeros(length(sensor_data_sgmntd{1}),1);

    sensor_data_sgmntd_atleast_1_type{1} = zeros(length(sensor_data_sgmntd{1}),1);
    sensor_data_sgmntd_atleast_2_type{1} = zeros(length(sensor_data_sgmntd{1}),1);
    sensor_data_sgmntd_atleast_3_type{1} = zeros(length(sensor_data_sgmntd{1}),1);

    if (n_label) % When there is a detection by the sensor system

        for k = 1 : n_label

            L_min = find(sensor_data_sgmntd_cmbd_all_sensors_labeled == k, 1 ); % Sample no. corresponding to start of the label
            L_max = find(sensor_data_sgmntd_cmbd_all_sensors_labeled == k, 1, 'last' ); % Sample no. corresponding to end of the label

            indv_detection_map = zeros(length(sensor_data_sgmntd{1}),1); % Need to be initialized before every detection matching
            indv_detection_map(L_min:L_max) = 1; % mapping individual sensation data

            % For detection by at least n sensors
            tmp_var = 0; % Variable to hold number of common sensors for each detection
            for j = 1:n_FM_sensors
                if (sum(indv_detection_map.*sensor_data_sgmntd{j})) % Non-zero value indicates intersection
                    tmp_var = tmp_var + 1;
                end
            end

            switch tmp_var
                case 1
                    sensor_data_sgmntd_atleast_1_sensor{1}(L_min:L_max) = 1;
                case 2
                    sensor_data_sgmntd_atleast_1_sensor{1}(L_min:L_max) = 1;
                    sensor_data_sgmntd_atleast_2_sensor{1}(L_min:L_max) = 1;
                case 3
                    sensor_data_sgmntd_atleast_1_sensor{1}(L_min:L_max) = 1;
                    sensor_data_sgmntd_atleast_2_sensor{1}(L_min:L_max) = 1;
                    sensor_data_sgmntd_atleast_3_sensor{1}(L_min:L_max) = 1;
                case 4
                    sensor_data_sgmntd_atleast_1_sensor{1}(L_min:L_max) = 1;
                    sensor_data_sgmntd_atleast_2_sensor{1}(L_min:L_max) = 1;
                    sensor_data_sgmntd_atleast_3_sensor{1}(L_min:L_max) = 1;
                    sensor_data_sgmntd_atleast_4_sensor{1}(L_min:L_max) = 1;
                case 5
                    sensor_data_sgmntd_atleast_1_sensor{1}(L_min:L_max) = 1;
                    sensor_data_sgmntd_atleast_2_sensor{1}(L_min:L_max) = 1;
                    sensor_data_sgmntd_atleast_3_sensor{1}(L_min:L_max) = 1;
                    sensor_data_sgmntd_atleast_4_sensor{1}(L_min:L_max) = 1;
                    sensor_data_sgmntd_atleast_5_sensor{1}(L_min:L_max) = 1;
                case 6
                    sensor_data_sgmntd_atleast_1_sensor{1}(L_min:L_max) = 1;
                    sensor_data_sgmntd_atleast_2_sensor{1}(L_min:L_max) = 1;
                    sensor_data_sgmntd_atleast_3_sensor{1}(L_min:L_max) = 1;
                    sensor_data_sgmntd_atleast_4_sensor{1}(L_min:L_max) = 1;
                    sensor_data_sgmntd_atleast_5_sensor{1}(L_min:L_max) = 1;
                    sensor_data_sgmntd_atleast_6_sensor{1}(L_min:L_max) = 1;
                otherwise
                    disp('This would never print.')
            end

            % For detection by at least n type of sensors
            tmp_var = 0; % Variable to hold number of common sensors for each detection
            for j = 1:n_FM_sensors/2
                if (sum(indv_detection_map.*sensor_data_sgmntd_cmbd_multi_type_sensors_OR{j})) % Non-zero value indicates intersection
                    tmp_var = tmp_var + 1;
                end
            end

            switch tmp_var
                case 1
                    sensor_data_sgmntd_atleast_1_type{1}(L_min:L_max) = 1;
                case 2
                    sensor_data_sgmntd_atleast_1_type{1}(L_min:L_max) = 1;
                    sensor_data_sgmntd_atleast_2_type{1}(L_min:L_max) = 1;
                case 3
                    sensor_data_sgmntd_atleast_1_type{1}(L_min:L_max) = 1;
                    sensor_data_sgmntd_atleast_2_type{1}(L_min:L_max) = 1;
                    sensor_data_sgmntd_atleast_3_type{1}(L_min:L_max) = 1;
                otherwise
                    disp('This would never print.')
            end
        end

    end
    %

    % ~~~~~~~~~~~~~~~~~ Matching with maternal sensation ~~~~~~~~~~~~~~~~~%
    %   match_with_m_sensation() function will be used here-
    %   Input variables:  sensor_data_sgmntd- A cell variable with single cell/multiple cells.
    %                                         Each cell contains data from a sensor or a combination.
    %                     sensation_data, IMU_map, M_sntn_Map- cell variables with single cell
    %                     ext_bakward, ext_forward, FM_dilation_time- scalar values
    %                     Fs_sensor, Fs_sensation- scalar values
    %   Output variables: TPD, FPD, TND, FND- vectors with number of rows equal to the
    %                                         number of cells in the sensor_data_sgmntd
    
    % For individual sensors
    %   The input argument 'sensor_data_sgmntd' will be a multi-cell 
    %   variable with number of cell = number of sensor data.
    %   Hence, the function will return 6 x 1 vectors for each output 
    %   argument because the input data has 6 cells.
    %   The values are then stored in cell variables to make it compatable
    %   with the performance analysis section.
    [TPD_indv_sensors, FPD_indv_sensors, TND_indv_sensors, FND_indv_sensors] = match_with_m_sensation(sensor_data_sgmntd, sensation_data_SD1_trimd{i}, ...
        IMU_map{i}, M_sntn_map{i}, ext_backward, ext_forward, FM_dilation_time, Fs_sensor, Fs_sensation);
    for j = 1 : n_FM_sensors
        TPD_all_indv_sensors_indv{j}(i,1) = TPD_indv_sensors(j); 
        FPD_all_indv_sensors_indv{j}(i,1) = FPD_indv_sensors(j);
        TND_all_indv_sensors_indv{j}(i,1) = TND_indv_sensors(j);
        FND_all_indv_sensors_indv{j}(i,1) = FND_indv_sensors(j);
        % Above variables will have 6 cell arrays, each for one sensors.
        % Each cell will contain a vector with row numbers = n_data_files.
    end
    %

    % For combined data
    %   input argument 'sensor_data_sgmntd' will be a single-cell variable. 
    %   Hence the function will return a scalar value for each output argument.
    %   Cell variable is used to store the output data to make it 
    %   compatible with the performance analysis section.
    %
    % Combined left and right sensors
    [TPD_Acstc_indv{1}(i,1), FPD_Acstc_indv{1}(i,1), TND_Acstc_indv{1}(i,1), FND_Acstc_indv{1}(i,1)] ...
        = match_with_m_sensation(sensor_data_sgmntd_Acstc, sensation_data_SD1_trimd{i}, IMU_map{i}, M_sntn_map{i}, ext_backward, ...
        ext_forward, FM_dilation_time, Fs_sensor, Fs_sensation); 
    [TPD_Aclm_indv{1}(i,1), FPD_Aclm_indv{1}(i,1), TND_Aclm_indv{1}(i,1), FND_Aclm_indv{1}(i,1)] ...
        = match_with_m_sensation(sensor_data_sgmntd_Aclm, sensation_data_SD1_trimd{i}, IMU_map{i}, M_sntn_map{i}, ext_backward, ...
        ext_forward, FM_dilation_time, Fs_sensor, Fs_sensation);
    [TPD_Pzplt_indv{1}(i,1), FPD_Pzplt_indv{1}(i,1), TND_Pzplt_indv{1}(i,1), FND_Pzplt_indv{1}(i,1)] ...
        = match_with_m_sensation(sensor_data_sgmntd_Pzplt, sensation_data_SD1_trimd{i}, IMU_map{i}, M_sntn_map{i}, ext_backward, ...
        ext_forward, FM_dilation_time, Fs_sensor, Fs_sensation);

    % Combination with sensor fusion scheme 1
    [TPD_Aclm_OR_Acstc_indv{1}(i,1), FPD_Aclm_OR_Acstc_indv{1}(i,1), TND_Aclm_OR_Acstc_indv{1}(i,1), FND_Aclm_OR_Acstc_indv{1}(i,1)] ...
        = match_with_m_sensation(sensor_data_sgmntd_Aclm_OR_Acstc, sensation_data_SD1_trimd{i}, IMU_map{i}, M_sntn_map{i}, ext_backward, ...
        ext_forward, FM_dilation_time, Fs_sensor, Fs_sensation);
    [TPD_Aclm_OR_Pzplt_indv{1}(i,1), FPD_Aclm_OR_Pzplt_indv{1}(i,1), TND_Aclm_OR_Pzplt_indv{1}(i,1), FND_Aclm_OR_Pzplt_indv{1}(i,1)] ...
        = match_with_m_sensation(sensor_data_sgmntd_Aclm_OR_Pzplt, sensation_data_SD1_trimd{i}, IMU_map{i}, M_sntn_map{i}, ext_backward, ...
        ext_forward, FM_dilation_time, Fs_sensor, Fs_sensation);
    [TPD_Acstc_OR_Pzplt_indv{1}(i,1), FPD_Acstc_OR_Pzplt_indv{1}(i,1), TND_Acstc_OR_Pzplt_indv{1}(i,1), FND_Acstc_OR_Pzplt_indv{1}(i,1)] ...
        = match_with_m_sensation(sensor_data_sgmntd_Acstc_OR_Pzplt, sensation_data_SD1_trimd{i}, IMU_map{i}, M_sntn_map{i}, ext_backward, ...
        ext_forward, FM_dilation_time, Fs_sensor, Fs_sensation);
    [TPD_all_sensors_OR_cmbd_indv{1}(i,1), FPD_all_sensors_OR_cmbd_indv{1}(i,1), TND_all_sensors_OR_cmbd_indv{1}(i,1), FND_all_sensors_OR_cmbd_indv{1}(i,1)] ...
        = match_with_m_sensation(sensor_data_sgmntd_cmbd_all_sensors_OR, sensation_data_SD1_trimd{i}, IMU_map{i}, M_sntn_map{i}, ext_backward, ...
        ext_forward, FM_dilation_time, Fs_sensor, Fs_sensation);
    %
    % Combination with sensor fusion scheme 2
    [TPD_Aclm_AND_Acstc_indv{1}(i,1), FPD_Aclm_AND_Acstc_indv{1}(i,1), TND_Aclm_AND_Acstc_indv{1}(i,1), FND_Aclm_AND_Acstc_indv{1}(i,1)] ...
        = match_with_m_sensation(sensor_data_sgmntd_Aclm_AND_Acstc, sensation_data_SD1_trimd{i}, IMU_map{i}, M_sntn_map{i}, ext_backward, ...
        ext_forward, FM_dilation_time, Fs_sensor, Fs_sensation);
    [TPD_Aclm_AND_Pzplt_indv{1}(i,1), FPD_Aclm_AND_Pzplt_indv{1}(i,1), TND_Aclm_AND_Pzplt_indv{1}(i,1), FND_Aclm_AND_Pzplt_indv{1}(i,1)] ...
        = match_with_m_sensation(sensor_data_sgmntd_Aclm_AND_Pzplt, sensation_data_SD1_trimd{i}, IMU_map{i}, M_sntn_map{i}, ext_backward, ...
        ext_forward, FM_dilation_time, Fs_sensor, Fs_sensation);
    [TPD_Acstc_AND_Pzplt_indv{1}(i,1), FPD_Acstc_AND_Pzplt_indv{1}(i,1), TND_Acstc_AND_Pzplt_indv{1}(i,1), FND_Acstc_AND_Pzplt_indv{1}(i,1)] ...
        = match_with_m_sensation(sensor_data_sgmntd_Acstc_AND_Pzplt, sensation_data_SD1_trimd{i}, IMU_map{i}, M_sntn_map{i}, ext_backward, ...
        ext_forward, FM_dilation_time, Fs_sensor, Fs_sensation);
    [TPD_all_sensors_AND_cmbd_indv{1}(i,1), FPD_all_sensors_AND_cmbd_indv{1}(i,1), TND_all_sensors_AND_cmbd_indv{1}(i,1), FND_all_sensors_AND_cmbd_indv{1}(i,1)] ...
        = match_with_m_sensation(sensor_data_sgmntd_cmbd_all_sensors_AND, sensation_data_SD1_trimd{i}, IMU_map{i}, M_sntn_map{i}, ext_backward, ...
        ext_forward, FM_dilation_time, Fs_sensor, Fs_sensation); 
    %
    % Combination with sensor fusion scheme 3
    [TPD_atleast_1_sensor_indv{1}(i,1), FPD_atleast_1_sensor_indv{1}(i,1), TND_atleast_1_sensor_indv{1}(i,1), FND_atleast_1_sensor_indv{1}(i,1)] ...
        = match_with_m_sensation(sensor_data_sgmntd_atleast_1_sensor, sensation_data_SD1_trimd{i}, IMU_map{i}, M_sntn_map{i}, ext_backward, ...
        ext_forward, FM_dilation_time, Fs_sensor, Fs_sensation); 
    [TPD_atleast_2_sensor_indv{1}(i,1), FPD_atleast_2_sensor_indv{1}(i,1), TND_atleast_2_sensor_indv{1}(i,1), FND_atleast_2_sensor_indv{1}(i,1)] ...
        = match_with_m_sensation(sensor_data_sgmntd_atleast_2_sensor, sensation_data_SD1_trimd{i}, IMU_map{i}, M_sntn_map{i}, ext_backward, ...
        ext_forward, FM_dilation_time, Fs_sensor, Fs_sensation);
    [TPD_atleast_3_sensor_indv{1}(i,1), FPD_atleast_3_sensor_indv{1}(i,1), TND_atleast_3_sensor_indv{1}(i,1), FND_atleast_3_sensor_indv{1}(i,1)] ...
        = match_with_m_sensation(sensor_data_sgmntd_atleast_3_sensor, sensation_data_SD1_trimd{i}, IMU_map{i}, M_sntn_map{i}, ext_backward, ...
        ext_forward, FM_dilation_time, Fs_sensor, Fs_sensation); 
    [TPD_atleast_4_sensor_indv{1}(i,1), FPD_atleast_4_sensor_indv{1}(i,1), TND_atleast_4_sensor_indv{1}(i,1), FND_atleast_4_sensor_indv{1}(i,1)] ...
        = match_with_m_sensation(sensor_data_sgmntd_atleast_4_sensor, sensation_data_SD1_trimd{i}, IMU_map{i}, M_sntn_map{i}, ext_backward, ...
        ext_forward, FM_dilation_time, Fs_sensor, Fs_sensation); 
    [TPD_atleast_5_sensor_indv{1}(i,1), FPD_atleast_5_sensor_indv{1}(i,1), TND_atleast_5_sensor_indv{1}(i,1), FND_atleast_5_sensor_indv{1}(i,1)] ...
        = match_with_m_sensation(sensor_data_sgmntd_atleast_5_sensor, sensation_data_SD1_trimd{i}, IMU_map{i}, M_sntn_map{i}, ext_backward, ...
        ext_forward, FM_dilation_time, Fs_sensor, Fs_sensation); 
    [TPD_atleast_6_sensor_indv{1}(i,1), FPD_atleast_6_sensor_indv{1}(i,1), TND_atleast_6_sensor_indv{1}(i,1), FND_atleast_6_sensor_indv{1}(i,1)] ...
        = match_with_m_sensation(sensor_data_sgmntd_atleast_6_sensor, sensation_data_SD1_trimd{i}, IMU_map{i}, M_sntn_map{i}, ext_backward, ...
        ext_forward, FM_dilation_time, Fs_sensor, Fs_sensation); 

    [TPD_atleast_1_type_indv{1}(i,1), FPD_atleast_1_type_indv{1}(i,1), TND_atleast_1_type_indv{1}(i,1), FND_atleast_1_type_indv{1}(i,1)] ...
        = match_with_m_sensation(sensor_data_sgmntd_atleast_1_type, sensation_data_SD1_trimd{i}, IMU_map{i}, M_sntn_map{i}, ext_backward, ...
        ext_forward, FM_dilation_time, Fs_sensor, Fs_sensation); 
    [TPD_atleast_2_type_indv{1}(i,1), FPD_atleast_2_type_indv{1}(i,1), TND_atleast_2_type_indv{1}(i,1), FND_atleast_2_type_indv{1}(i,1)] ...
        = match_with_m_sensation(sensor_data_sgmntd_atleast_2_type, sensation_data_SD1_trimd{i}, IMU_map{i}, M_sntn_map{i}, ext_backward, ...
        ext_forward, FM_dilation_time, Fs_sensor, Fs_sensation);
    [TPD_atleast_3_type_indv{1}(i,1), FPD_atleast_3_type_indv{1}(i,1), TND_atleast_3_type_indv{1}(i,1), FND_atleast_3_type_indv{1}(i,1)] ...
        = match_with_m_sensation(sensor_data_sgmntd_atleast_3_type, sensation_data_SD1_trimd{i}, IMU_map{i}, M_sntn_map{i}, ext_backward, ...
        ext_forward, FM_dilation_time, Fs_sensor, Fs_sensation);
    %
end
fprintf('Performance analysis completed.\n\n')
%

% ----------------------- Performance analysis ---------------------------%
% This section will use get_performance_params() function
%   Input variables:  TPD_all, FPD_all, TND_all, FND_all- single cell/multi-cell variable.
%                     Number of cell indicates number of sensor data or
%                     combination data provided together.
%                     Each cell containes a vector with row size = n_data_files 
%
%   Output variables: SEN_all,PPV_all,SPE_all,ACC_all,FS_all,FPR_all- cell variable with
%                     size same as the input variables.   

% For individual data sets                  
[SEN_all_indv_sensors_indv, PPV_all_indv_sensors_indv, SPE_all_indv_sensors_indv, ACC_all_indv_sensors_indv, FS_all_indv_sensors_indv, FPR_all_indv_sensors_indv] ...
    = get_performance_params(TPD_all_indv_sensors_indv, FPD_all_indv_sensors_indv, TND_all_indv_sensors_indv, FND_all_indv_sensors_indv);

[SEN_Aclm_indv, PPV_Aclm_indv, SPE_Aclm_indv, ACC_Aclm_indv, FS_Aclm_indv, FPR_Aclm_indv] ...
    = get_performance_params(TPD_Aclm_indv, FPD_Aclm_indv, TND_Aclm_indv, FND_Aclm_indv);
[SEN_Acstc_indv, PPV_Acstc_indv, SPE_Acstc_indv, ACC_Acstc_indv, FS_Acstc_indv, FPR_Acstc_indv] ...
    = get_performance_params(TPD_Acstc_indv, FPD_Acstc_indv, TND_Acstc_indv, FND_Acstc_indv); 
[SEN_Pzplt_indv, PPV_Pzplt_indv, SPE_Pzplt_indv, ACC_Pzplt_indv, FS_Pzplt_indv, FPR_Pzplt_indv] ...
    = get_performance_params(TPD_Pzplt_indv, FPD_Pzplt_indv, TND_Pzplt_indv, FND_Pzplt_indv);

[SEN_Aclm_OR_Acstc_indv, PPV_Aclm_OR_Acstc_indv, SPE_Aclm_OR_Acstc_indv, ACC_Aclm_OR_Acstc_indv, FS_Aclm_OR_Acstc_indv, FPR_Aclm_OR_Acstc_indv] ...
    = get_performance_params(TPD_Aclm_OR_Acstc_indv, FPD_Aclm_OR_Acstc_indv, TND_Aclm_OR_Acstc_indv, FND_Aclm_OR_Acstc_indv);
[SEN_Aclm_OR_Pzplt_indv, PPV_Aclm_OR_Pzplt_indv, SPE_Aclm_OR_Pzplt_indv, ACC_Aclm_OR_Pzplt_indv, FS_Aclm_OR_Pzplt_indv, FPR_Aclm_OR_Pzplt_indv] ...
    = get_performance_params(TPD_Aclm_OR_Pzplt_indv, FPD_Aclm_OR_Pzplt_indv, TND_Aclm_OR_Pzplt_indv, FND_Aclm_OR_Pzplt_indv);
[SEN_Acstc_OR_Pzplt_indv, PPV_Acstc_OR_Pzplt_indv, SPE_Acstc_OR_Pzplt_indv, ACC_Acstc_OR_Pzplt_indv, FS_Acstc_OR_Pzplt_indv, FPR_Acstc_OR_Pzplt_indv] ...
    = get_performance_params(TPD_Acstc_OR_Pzplt_indv, FPD_Acstc_OR_Pzplt_indv, TND_Acstc_OR_Pzplt_indv, FND_Acstc_OR_Pzplt_indv);
[SEN_all_sensors_OR_cmbd_indv, PPV_all_sensors_OR_cmbd_indv, SPE_all_sensors_OR_cmbd_indv, ACC_all_sensors_OR_cmbd_indv, FS_all_sensors_OR_cmbd_indv, FPR_all_OR_sensors_cmbd_indv] ...
    =  get_performance_params(TPD_all_sensors_OR_cmbd_indv, FPD_all_sensors_OR_cmbd_indv, TND_all_sensors_OR_cmbd_indv, FND_all_sensors_OR_cmbd_indv);

[SEN_Aclm_AND_Acstc_indv, PPV_Aclm_AND_Acstc_indv, SPE_Aclm_AND_Acstc_indv, ACC_Aclm_AND_Acstc_indv, FS_Aclm_AND_Acstc_indv, FPR_Aclm_AND_Acstc_indv] ...
    = get_performance_params(TPD_Aclm_AND_Acstc_indv, FPD_Aclm_AND_Acstc_indv, TND_Aclm_AND_Acstc_indv, FND_Aclm_AND_Acstc_indv);
[SEN_Aclm_AND_Pzplt_indv, PPV_Aclm_AND_Pzplt_indv, SPE_Aclm_AND_Pzplt_indv, ACC_Aclm_AND_Pzplt_indv, FS_Aclm_AND_Pzplt_indv, FPR_Aclm_AND_Pzplt_indv] ...
    = get_performance_params(TPD_Aclm_AND_Pzplt_indv, FPD_Aclm_AND_Pzplt_indv, TND_Aclm_AND_Pzplt_indv, FND_Aclm_AND_Pzplt_indv);
[SEN_Acstc_AND_Pzplt_indv, PPV_Acstc_AND_Pzplt_indv, SPE_Acstc_AND_Pzplt_indv, ACC_Acstc_AND_Pzplt_indv, FS_Acstc_AND_Pzplt_indv, FPR_Acstc_AND_Pzplt_indv] ...
    = get_performance_params(TPD_Acstc_AND_Pzplt_indv, FPD_Acstc_AND_Pzplt_indv, TND_Acstc_AND_Pzplt_indv, FND_Acstc_AND_Pzplt_indv);
[SEN_all_sensors_AND_cmbd_indv, PPV_all_sensors_AND_cmbd_indv, SPE_all_sensors_AND_cmbd_indv, ACC_all_sensors_AND_cmbd_indv, FS_all_sensors_AND_cmbd_indv, FPR_all_AND_sensors_cmbd_indv] ...
    =  get_performance_params(TPD_all_sensors_AND_cmbd_indv, FPD_all_sensors_AND_cmbd_indv, TND_all_sensors_AND_cmbd_indv, FND_all_sensors_AND_cmbd_indv);

[SEN_atleast_1_sensor_indv, PPV_atleast_1_sensor_indv, SPE_atleast_1_sensor_indv, ACC_atleast_1_sensor_indv, FS_atleast_1_sensor_indv, FPR_atleast_1_sensor_indv] ...
    = get_performance_params(TPD_atleast_1_sensor_indv, FPD_atleast_1_sensor_indv, TND_atleast_1_sensor_indv, FND_atleast_1_sensor_indv);
[SEN_atleast_2_sensor_indv, PPV_atleast_2_sensor_indv, SPE_atleast_2_sensor_indv, ACC_atleast_2_sensor_indv, FS_atleast_2_sensor_indv, FPR_atleast_2_sensor_indv] ...
    = get_performance_params(TPD_atleast_2_sensor_indv, FPD_atleast_2_sensor_indv, TND_atleast_2_sensor_indv, FND_atleast_2_sensor_indv);
[SEN_atleast_3_sensor_indv, PPV_atleast_3_sensor_indv, SPE_atleast_3_sensor_indv, ACC_atleast_3_sensor_indv, FS_atleast_3_sensor_indv, FPR_atleast_3_sensor_indv] ...
    = get_performance_params(TPD_atleast_3_sensor_indv, FPD_atleast_3_sensor_indv, TND_atleast_3_sensor_indv, FND_atleast_3_sensor_indv);
[SEN_atleast_4_sensor_indv, PPV_atleast_4_sensor_indv, SPE_atleast_4_sensor_indv, ACC_atleast_4_sensor_indv, FS_atleast_4_sensor_indv, FPR_atleast_4_sensor_indv] ...
    = get_performance_params(TPD_atleast_4_sensor_indv, FPD_atleast_4_sensor_indv, TND_atleast_4_sensor_indv, FND_atleast_4_sensor_indv);
[SEN_atleast_5_sensor_indv, PPV_atleast_5_sensor_indv, SPE_atleast_5_sensor_indv, ACC_atleast_5_sensor_indv, FS_atleast_5_sensor_indv, FPR_atleast_5_sensor_indv] ...
    = get_performance_params(TPD_atleast_5_sensor_indv, FPD_atleast_5_sensor_indv, TND_atleast_5_sensor_indv, FND_atleast_5_sensor_indv);
[SEN_atleast_6_sensor_indv, PPV_atleast_6_sensor_indv, SPE_atleast_6_sensor_indv, ACC_atleast_6_sensor_indv, FS_atleast_6_sensor_indv, FPR_atleast_6_sensor_indv] ...
    = get_performance_params(TPD_atleast_6_sensor_indv, FPD_atleast_6_sensor_indv, TND_atleast_6_sensor_indv, FND_atleast_6_sensor_indv);

[SEN_atleast_1_type_indv, PPV_atleast_1_type_indv, SPE_atleast_1_type_indv, ACC_atleast_1_type_indv, FS_atleast_1_type_indv, FPR_atleast_1_type_indv] ...
    = get_performance_params(TPD_atleast_1_type_indv, FPD_atleast_1_type_indv, TND_atleast_1_type_indv, FND_atleast_1_type_indv);
[SEN_atleast_2_type_indv, PPV_atleast_2_type_indv, SPE_atleast_2_type_indv, ACC_atleast_2_type_indv, FS_atleast_2_type_indv, FPR_atleast_2_type_indv] ...
    = get_performance_params(TPD_atleast_2_type_indv, FPD_atleast_2_type_indv, TND_atleast_2_type_indv, FND_atleast_2_type_indv);
[SEN_atleast_3_type_indv, PPV_atleast_3_type_indv, SPE_atleast_3_type_indv, ACC_atleast_3_type_indv, FS_atleast_3_type_indv, FPR_atleast_3_type_indv] ...
    = get_performance_params(TPD_atleast_3_type_indv, FPD_atleast_3_type_indv, TND_atleast_3_type_indv, FND_atleast_3_type_indv);

% For all the data sets combinedly
TPD_all_indv_sensors_overall = cell(1, n_FM_sensors);
FPD_all_indv_sensors_overall = cell(1, n_FM_sensors);
TND_all_indv_sensors_overall = cell(1, n_FM_sensors);
FND_all_indv_sensors_overall = cell(1, n_FM_sensors);

for j = 1 : n_FM_sensors    
    TPD_all_indv_sensors_overall{j} = sum(TPD_all_indv_sensors_indv{j} ,1); % Sums up the elements of each column and returns a row vector
    FPD_all_indv_sensors_overall{j} = sum(FPD_all_indv_sensors_indv{j} ,1);
    TND_all_indv_sensors_overall{j} = sum(TND_all_indv_sensors_indv{j} ,1);
    FND_all_indv_sensors_overall{j} = sum(FND_all_indv_sensors_indv{j} ,1);
end

TPD_Aclm_overall{1}  = sum(TPD_Aclm_indv{1} ,1); % Sums up across a column 
FPD_Aclm_overall{1}  = sum(FPD_Aclm_indv{1} ,1);
TND_Aclm_overall{1}  = sum(TND_Aclm_indv{1} ,1);
FND_Aclm_overall{1}  = sum(FND_Aclm_indv{1} ,1);
TPD_Acstc_overall{1} = sum(TPD_Acstc_indv{1} ,1); 
FPD_Acstc_overall{1} = sum(FPD_Acstc_indv{1} ,1);
TND_Acstc_overall{1} = sum(TND_Acstc_indv{1} ,1);
FND_Acstc_overall{1} = sum(FND_Acstc_indv{1} ,1);
TPD_Pzplt_overall{1} = sum(TPD_Pzplt_indv{1} ,1);
FPD_Pzplt_overall{1} = sum(FPD_Pzplt_indv{1} ,1);
TND_Pzplt_overall{1} = sum(TND_Pzplt_indv{1} ,1);
FND_Pzplt_overall{1} = sum(FND_Pzplt_indv{1} ,1);

TPD_Aclm_OR_Acstc_overall{1}  = sum(TPD_Aclm_OR_Acstc_indv{1} ,1);
FPD_Aclm_OR_Acstc_overall{1}  = sum(FPD_Aclm_OR_Acstc_indv{1} ,1);
TND_Aclm_OR_Acstc_overall{1}  = sum(TND_Aclm_OR_Acstc_indv{1} ,1);
FND_Aclm_OR_Acstc_overall{1}  = sum(FND_Aclm_OR_Acstc_indv{1} ,1);
TPD_Aclm_OR_Pzplt_overall{1}  = sum(TPD_Aclm_OR_Pzplt_indv{1} ,1);
FPD_Aclm_OR_Pzplt_overall{1}  = sum(FPD_Aclm_OR_Pzplt_indv{1} ,1);
TND_Aclm_OR_Pzplt_overall{1}  = sum(TND_Aclm_OR_Pzplt_indv{1} ,1);
FND_Aclm_OR_Pzplt_overall{1}  = sum(FND_Aclm_OR_Pzplt_indv{1} ,1); 
TPD_Acstc_OR_Pzplt_overall{1} = sum(TPD_Acstc_OR_Pzplt_indv{1} ,1);
FPD_Acstc_OR_Pzplt_overall{1} = sum(FPD_Acstc_OR_Pzplt_indv{1} ,1);
TND_Acstc_OR_Pzplt_overall{1} = sum(TND_Acstc_OR_Pzplt_indv{1} ,1);
FND_Acstc_OR_Pzplt_overall{1} = sum(FND_Acstc_OR_Pzplt_indv{1} ,1);
TPD_all_sensors_OR_cmbd_overall{1} = sum(TPD_all_sensors_OR_cmbd_indv{1} ,1);
FPD_all_sensors_OR_cmbd_overall{1} = sum(FPD_all_sensors_OR_cmbd_indv{1} ,1);
TND_all_sensors_OR_cmbd_overall{1} = sum(TND_all_sensors_OR_cmbd_indv{1} ,1);
FND_all_sensors_OR_cmbd_overall{1} = sum(FND_all_sensors_OR_cmbd_indv{1} ,1);

TPD_Aclm_AND_Acstc_overall{1}  = sum(TPD_Aclm_AND_Acstc_indv{1} ,1);
FPD_Aclm_AND_Acstc_overall{1}  = sum(FPD_Aclm_AND_Acstc_indv{1} ,1);
TND_Aclm_AND_Acstc_overall{1}  = sum(TND_Aclm_AND_Acstc_indv{1} ,1);
FND_Aclm_AND_Acstc_overall{1}  = sum(FND_Aclm_AND_Acstc_indv{1} ,1);
TPD_Aclm_AND_Pzplt_overall{1}  = sum(TPD_Aclm_AND_Pzplt_indv{1} ,1);
FPD_Aclm_AND_Pzplt_overall{1}  = sum(FPD_Aclm_AND_Pzplt_indv{1} ,1);
TND_Aclm_AND_Pzplt_overall{1}  = sum(TND_Aclm_AND_Pzplt_indv{1} ,1);
FND_Aclm_AND_Pzplt_overall{1}  = sum(FND_Aclm_AND_Pzplt_indv{1} ,1); 
TPD_Acstc_AND_Pzplt_overall{1} = sum(TPD_Acstc_AND_Pzplt_indv{1} ,1);
FPD_Acstc_AND_Pzplt_overall{1} = sum(FPD_Acstc_AND_Pzplt_indv{1} ,1);
TND_Acstc_AND_Pzplt_overall{1} = sum(TND_Acstc_AND_Pzplt_indv{1} ,1);
FND_Acstc_AND_Pzplt_overall{1} = sum(FND_Acstc_AND_Pzplt_indv{1} ,1);
TPD_all_sensors_AND_cmbd_overall{1} = sum(TPD_all_sensors_AND_cmbd_indv{1} ,1);
FPD_all_sensors_AND_cmbd_overall{1} = sum(FPD_all_sensors_AND_cmbd_indv{1} ,1);
TND_all_sensors_AND_cmbd_overall{1} = sum(TND_all_sensors_AND_cmbd_indv{1} ,1);
FND_all_sensors_AND_cmbd_overall{1} = sum(FND_all_sensors_AND_cmbd_indv{1} ,1);

TPD_atleast_1_sensor_overall{1} = sum(TPD_atleast_1_sensor_indv{1} ,1); % Sums up across a column 
FPD_atleast_1_sensor_overall{1} = sum(FPD_atleast_1_sensor_indv{1} ,1);
TND_atleast_1_sensor_overall{1} = sum(TND_atleast_1_sensor_indv{1} ,1);
FND_atleast_1_sensor_overall{1} = sum(FND_atleast_1_sensor_indv{1} ,1);
TPD_atleast_2_sensor_overall{1} = sum(TPD_atleast_2_sensor_indv{1} ,1); 
FPD_atleast_2_sensor_overall{1} = sum(FPD_atleast_2_sensor_indv{1} ,1);
TND_atleast_2_sensor_overall{1} = sum(TND_atleast_2_sensor_indv{1} ,1);
FND_atleast_2_sensor_overall{1} = sum(FND_atleast_2_sensor_indv{1} ,1);
TPD_atleast_3_sensor_overall{1} = sum(TPD_atleast_3_sensor_indv{1} ,1); 
FPD_atleast_3_sensor_overall{1} = sum(FPD_atleast_3_sensor_indv{1} ,1);
TND_atleast_3_sensor_overall{1} = sum(TND_atleast_3_sensor_indv{1} ,1);
FND_atleast_3_sensor_overall{1} = sum(FND_atleast_3_sensor_indv{1} ,1);
TPD_atleast_4_sensor_overall{1} = sum(TPD_atleast_4_sensor_indv{1} ,1); 
FPD_atleast_4_sensor_overall{1} = sum(FPD_atleast_4_sensor_indv{1} ,1);
TND_atleast_4_sensor_overall{1} = sum(TND_atleast_4_sensor_indv{1} ,1);
FND_atleast_4_sensor_overall{1} = sum(FND_atleast_4_sensor_indv{1} ,1);
TPD_atleast_5_sensor_overall{1} = sum(TPD_atleast_5_sensor_indv{1} ,1); 
FPD_atleast_5_sensor_overall{1} = sum(FPD_atleast_5_sensor_indv{1} ,1);
TND_atleast_5_sensor_overall{1} = sum(TND_atleast_5_sensor_indv{1} ,1);
FND_atleast_5_sensor_overall{1} = sum(FND_atleast_5_sensor_indv{1} ,1);
TPD_atleast_6_sensor_overall{1} = sum(TPD_atleast_6_sensor_indv{1} ,1); 
FPD_atleast_6_sensor_overall{1} = sum(FPD_atleast_6_sensor_indv{1} ,1);
TND_atleast_6_sensor_overall{1} = sum(TND_atleast_6_sensor_indv{1} ,1);
FND_atleast_6_sensor_overall{1} = sum(FND_atleast_6_sensor_indv{1} ,1);

TPD_atleast_1_type_overall{1} = sum(TPD_atleast_1_type_indv{1} ,1); % Sums up across a column 
FPD_atleast_1_type_overall{1} = sum(FPD_atleast_1_type_indv{1} ,1);
TND_atleast_1_type_overall{1} = sum(TND_atleast_1_type_indv{1} ,1);
FND_atleast_1_type_overall{1} = sum(FND_atleast_1_type_indv{1} ,1);
TPD_atleast_2_type_overall{1} = sum(TPD_atleast_2_type_indv{1} ,1); 
FPD_atleast_2_type_overall{1} = sum(FPD_atleast_2_type_indv{1} ,1);
TND_atleast_2_type_overall{1} = sum(TND_atleast_2_type_indv{1} ,1);
FND_atleast_2_type_overall{1} = sum(FND_atleast_2_type_indv{1} ,1);
TPD_atleast_3_type_overall{1} = sum(TPD_atleast_3_type_indv{1} ,1); 
FPD_atleast_3_type_overall{1} = sum(FPD_atleast_3_type_indv{1} ,1);
TND_atleast_3_type_overall{1} = sum(TND_atleast_3_type_indv{1} ,1);
FND_atleast_3_type_overall{1} = sum(FND_atleast_3_type_indv{1} ,1);

[SEN_all_indv_sensors_overall, PPV_all_indv_sensors_overall, SPE_all_indv_sensors_overall, ACC_all_indv_sensors_overall, FS_all_indv_sensors_overall, FPR_all_indv_sensors_overall] ...
    = get_performance_params(TPD_all_indv_sensors_overall, FPD_all_indv_sensors_overall, TND_all_indv_sensors_overall, FND_all_indv_sensors_overall); % Values for each sensor will be stored in each cell

[SEN_Acstc_overall, PPV_Acstc_overall, SPE_Acstc_overall, ACC_Acstc_overall, FS_Acstc_overall, FPR_Acstc_overall] ...
    = get_performance_params(TPD_Acstc_overall, FPD_Acstc_overall, TND_Acstc_overall, FND_Acstc_overall);
[SEN_Aclm_overall, PPV_Aclm_overall, SPE_Aclm_overall, ACC_Aclm_overall, FS_Aclm_overall, FPR_Aclm_overall] ...
    = get_performance_params(TPD_Aclm_overall, FPD_Aclm_overall, TND_Aclm_overall, FND_Aclm_overall);
[SEN_Pzplt_overall, PPV_Pzplt_overall, SPE_Pzplt_overall, ACC_Pzplt_overall, FS_Pzplt_overall, FPR_Pzplt_overall] ...
    = get_performance_params(TPD_Pzplt_overall, FPD_Pzplt_overall, TND_Pzplt_overall, FND_Pzplt_overall);
[SEN_Aclm_OR_Acstc_overall, PPV_Aclm_OR_Acstc_overall, SPE_Aclm_OR_Acstc_overall, ACC_Aclm_OR_Acstc_overall, FS_Aclm_OR_Acstc_overall, FPR_Aclm_OR_Acstc_overall] ...
    = get_performance_params(TPD_Aclm_OR_Acstc_overall, FPD_Aclm_OR_Acstc_overall, TND_Aclm_OR_Acstc_overall, FND_Aclm_OR_Acstc_overall);
[SEN_Aclm_OR_Pzplt_overall, PPV_Aclm_OR_Pzplt_overall, SPE_Aclm_OR_Pzplt_overall, ACC_Aclm_OR_Pzplt_overall, FS_Aclm_OR_Pzplt_overall, FPR_Aclm_OR_Pzplt_overall] ...
    = get_performance_params(TPD_Aclm_OR_Pzplt_overall, FPD_Aclm_OR_Pzplt_overall, TND_Aclm_OR_Pzplt_overall, FND_Aclm_OR_Pzplt_overall);
[SEN_Acstc_OR_Pzplt_overall, PPV_Acstc_OR_Pzplt_overall, SPE_Acstc_OR_Pzplt_overall, ACC_Acstc_OR_Pzplt_overall, FS_Acstc_OR_Pzplt_overall, FPR_Acstc_OR_Pzplt_overall] ...
    = get_performance_params(TPD_Acstc_OR_Pzplt_overall, FPD_Acstc_OR_Pzplt_overall, TND_Acstc_OR_Pzplt_overall, FND_Acstc_OR_Pzplt_overall);
[SEN_all_sensors_OR_cmbd_overall, PPV_all_sensors_OR_cmbd_overall, SPE_all_sensors_OR_cmbd_overall, ACC_all_sensors_OR_cmbd_overall, FS_all_sensors_OR_cmbd_overall, FPR_all_sensors_OR_cmbd_overall] ...
    = get_performance_params(TPD_all_sensors_OR_cmbd_overall, FPD_all_sensors_OR_cmbd_overall, TND_all_sensors_OR_cmbd_overall, FND_all_sensors_OR_cmbd_overall);

[SEN_Aclm_AND_Acstc_overall, PPV_Aclm_AND_Acstc_overall, SPE_Aclm_AND_Acstc_overall, ACC_Aclm_AND_Acstc_overall, FS_Aclm_AND_Acstc_overall, FPR_Aclm_AND_Acstc_overall] ...
    = get_performance_params(TPD_Aclm_AND_Acstc_overall, FPD_Aclm_AND_Acstc_overall, TND_Aclm_AND_Acstc_overall, FND_Aclm_AND_Acstc_overall);
[SEN_Aclm_AND_Pzplt_overall, PPV_Aclm_AND_Pzplt_overall, SPE_Aclm_AND_Pzplt_overall, ACC_Aclm_AND_Pzplt_overall, FS_Aclm_AND_Pzplt_overall, FPR_Aclm_AND_Pzplt_overall] ...
    = get_performance_params(TPD_Aclm_AND_Pzplt_overall, FPD_Aclm_AND_Pzplt_overall, TND_Aclm_AND_Pzplt_overall, FND_Aclm_AND_Pzplt_overall);
[SEN_Acstc_AND_Pzplt_overall, PPV_Acstc_AND_Pzplt_overall, SPE_Acstc_AND_Pzplt_overall, ACC_Acstc_AND_Pzplt_overall, FS_Acstc_AND_Pzplt_overall, FPR_Acstc_AND_Pzplt_overall] ...
    = get_performance_params(TPD_Acstc_AND_Pzplt_overall, FPD_Acstc_AND_Pzplt_overall, TND_Acstc_AND_Pzplt_overall, FND_Acstc_AND_Pzplt_overall);
[SEN_all_sensors_AND_cmbd_overall, PPV_all_sensors_AND_cmbd_overall, SPE_all_sensors_AND_cmbd_overall, ACC_all_sensors_AND_cmbd_overall, FS_all_sensors_AND_cmbd_overall, FPR_all_sensors_AND_cmbd_overall] ...
    = get_performance_params(TPD_all_sensors_AND_cmbd_overall, FPD_all_sensors_AND_cmbd_overall, TND_all_sensors_AND_cmbd_overall, FND_all_sensors_AND_cmbd_overall);

[SEN_atleast_1_sensor_overall, PPV_atleast_1_sensor_overall, SPE_atleast_1_sensor_overall, ACC_atleast_1_sensor_overall, FS_atleast_1_sensor_overall, FPR_atleast_1_sensor_overall] ...
    = get_performance_params(TPD_atleast_1_sensor_overall, FPD_atleast_1_sensor_overall, TND_atleast_1_sensor_overall, FND_atleast_1_sensor_overall);
[SEN_atleast_2_sensor_overall, PPV_atleast_2_sensor_overall, SPE_atleast_2_sensor_overall, ACC_atleast_2_sensor_overall, FS_atleast_2_sensor_overall, FPR_atleast_2_sensor_overall] ...
    = get_performance_params(TPD_atleast_2_sensor_overall, FPD_atleast_2_sensor_overall, TND_atleast_2_sensor_overall, FND_atleast_2_sensor_overall);
[SEN_atleast_3_sensor_overall, PPV_atleast_3_sensor_overall, SPE_atleast_3_sensor_overall, ACC_atleast_3_sensor_overall, FS_atleast_3_sensor_overall, FPR_atleast_3_sensor_overall] ...
    = get_performance_params(TPD_atleast_3_sensor_overall, FPD_atleast_3_sensor_overall, TND_atleast_3_sensor_overall, FND_atleast_3_sensor_overall);
[SEN_atleast_4_sensor_overall, PPV_atleast_4_sensor_overall, SPE_atleast_4_sensor_overall, ACC_atleast_4_sensor_overall, FS_atleast_4_sensor_overall, FPR_atleast_4_sensor_overall] ...
    = get_performance_params(TPD_atleast_4_sensor_overall, FPD_atleast_4_sensor_overall, TND_atleast_4_sensor_overall, FND_atleast_4_sensor_overall);
[SEN_atleast_5_sensor_overall, PPV_atleast_5_sensor_overall, SPE_atleast_5_sensor_overall, ACC_atleast_5_sensor_overall, FS_atleast_5_sensor_overall, FPR_atleast_5_sensor_overall] ...
    = get_performance_params(TPD_atleast_5_sensor_overall, FPD_atleast_5_sensor_overall, TND_atleast_5_sensor_overall, FND_atleast_5_sensor_overall);
[SEN_atleast_6_sensor_overall, PPV_atleast_6_sensor_overall, SPE_atleast_6_sensor_overall, ACC_atleast_6_sensor_overall, FS_atleast_6_sensor_overall, FPR_atleast_6_sensor_overall] ...
    = get_performance_params(TPD_atleast_6_sensor_overall, FPD_atleast_6_sensor_overall, TND_atleast_6_sensor_overall, FND_atleast_6_sensor_overall);

[SEN_atleast_1_type_overall, PPV_atleast_1_type_overall, SPE_atleast_1_type_overall, ACC_atleast_1_type_overall, FS_atleast_1_type_overall, FPR_atleast_1_type_overall] ...
    = get_performance_params(TPD_atleast_1_type_overall, FPD_atleast_1_type_overall, TND_atleast_1_type_overall, FND_atleast_1_type_overall);
[SEN_atleast_2_type_overall, PPV_atleast_2_type_overall, SPE_atleast_2_type_overall, ACC_atleast_2_type_overall, FS_atleast_2_type_overall, FPR_atleast_2_type_overall] ...
    = get_performance_params(TPD_atleast_2_type_overall, FPD_atleast_2_type_overall, TND_atleast_2_type_overall, FND_atleast_2_type_overall);
[SEN_atleast_3_type_overall, PPV_atleast_3_type_overall, SPE_atleast_3_type_overall, ACC_atleast_3_type_overall, FS_atleast_3_type_overall, FPR_atleast_3_type_overall] ...
    = get_performance_params(TPD_atleast_3_type_overall, FPD_atleast_3_type_overall, TND_atleast_3_type_overall, FND_atleast_3_type_overall);
%

% ====================== Storing data in a table =========================
% Determination of duration of each data files
duration_data_files = zeros(n_data_files, 1);

for i = 1 : n_data_files
    duration_data_files(i,1) = length(Acstc_data1{i})/(Fs_sensor*60); % Duration in minutes of the original data file before trimming
end
%

% Combining information for differen sensors
sensor_combinations = ["Left Accelerometer", "Right Accelerometer", "Left Acoustic", "Right Acoustic", "Left Piezo", "Right Piezo", "Combined Accelerometers", "Combined Acoustics", "Combined Piezos", ...
    "Accelerometers OR Acoustics",  "Accelerometers OR Piezoelectrics", "Acoustics OR Piezoelectrics",  "All sensors OR", ...
    "Accelerometers AND Acoustics", "Accelerometers AND Piezoelectrics", "Acoustics AND Piezoelectrics", "All sensors AND", ...
    "At least 1 sensor",         "At least 2 sensors",           "At least 3 sensors",   "At least 4 sensors",           "At least 5 sensors",           "At least 6 sensors", ...
    "At least 1 type of sensor", "At least 2 types of sensor",   "At least 3 types of sensor"];

FM_min_SN_used = "indv sensor type based";

FS_indv_total  = [FS_all_indv_sensors_indv, FS_Aclm_indv, FS_Acstc_indv, FS_Pzplt_indv,...
    FS_Aclm_OR_Acstc_indv,    FS_Aclm_OR_Pzplt_indv,    FS_Acstc_OR_Pzplt_indv,   FS_all_sensors_OR_cmbd_indv, ...
    FS_Aclm_AND_Acstc_indv,   FS_Aclm_AND_Pzplt_indv,   FS_Acstc_AND_Pzplt_indv,  FS_all_sensors_AND_cmbd_indv, ...
    FS_atleast_1_sensor_indv, FS_atleast_2_sensor_indv, FS_atleast_3_sensor_indv, FS_atleast_4_sensor_indv,     FS_atleast_5_sensor_indv, FS_atleast_6_sensor_indv, ...
    FS_atleast_1_type_indv,   FS_atleast_2_type_indv,   FS_atleast_3_type_indv];
SEN_indv_total  = [SEN_all_indv_sensors_indv, SEN_Aclm_indv, SEN_Acstc_indv, SEN_Pzplt_indv,...
    SEN_Aclm_OR_Acstc_indv,    SEN_Aclm_OR_Pzplt_indv,    SEN_Acstc_OR_Pzplt_indv,   SEN_all_sensors_OR_cmbd_indv, ...
    SEN_Aclm_AND_Acstc_indv,   SEN_Aclm_AND_Pzplt_indv,   SEN_Acstc_AND_Pzplt_indv,  SEN_all_sensors_AND_cmbd_indv, ...
    SEN_atleast_1_sensor_indv, SEN_atleast_2_sensor_indv, SEN_atleast_3_sensor_indv, SEN_atleast_4_sensor_indv,     SEN_atleast_5_sensor_indv, SEN_atleast_6_sensor_indv, ...
    SEN_atleast_1_type_indv,   SEN_atleast_2_type_indv,   SEN_atleast_3_type_indv];
PPV_indv_total  = [PPV_all_indv_sensors_indv, PPV_Aclm_indv, PPV_Acstc_indv, PPV_Pzplt_indv,...
    PPV_Aclm_OR_Acstc_indv,    PPV_Aclm_OR_Pzplt_indv,    PPV_Acstc_OR_Pzplt_indv,   PPV_all_sensors_OR_cmbd_indv, ...
    PPV_Aclm_AND_Acstc_indv,   PPV_Aclm_AND_Pzplt_indv,   PPV_Acstc_AND_Pzplt_indv,  PPV_all_sensors_AND_cmbd_indv, ...
    PPV_atleast_1_sensor_indv, PPV_atleast_2_sensor_indv, PPV_atleast_3_sensor_indv, PPV_atleast_4_sensor_indv,     PPV_atleast_5_sensor_indv, PPV_atleast_6_sensor_indv, ...
    PPV_atleast_1_type_indv,   PPV_atleast_2_type_indv,   PPV_atleast_3_type_indv];
SPE_indv_total  = [SPE_all_indv_sensors_indv, SPE_Aclm_indv, SPE_Acstc_indv, SPE_Pzplt_indv,...
    SPE_Aclm_OR_Acstc_indv,    SPE_Aclm_OR_Pzplt_indv,    SPE_Acstc_OR_Pzplt_indv,   SPE_all_sensors_OR_cmbd_indv, ...
    SPE_Aclm_AND_Acstc_indv,   SPE_Aclm_AND_Pzplt_indv,   SPE_Acstc_AND_Pzplt_indv,  SPE_all_sensors_AND_cmbd_indv, ...
    SPE_atleast_1_sensor_indv, SPE_atleast_2_sensor_indv, SPE_atleast_3_sensor_indv, SPE_atleast_4_sensor_indv,     SPE_atleast_5_sensor_indv, SPE_atleast_6_sensor_indv, ...
    SPE_atleast_1_type_indv,   SPE_atleast_2_type_indv,   SPE_atleast_3_type_indv];
ACC_indv_total  = [ACC_all_indv_sensors_indv, ACC_Aclm_indv, ACC_Acstc_indv, ACC_Pzplt_indv,...
    ACC_Aclm_OR_Acstc_indv,    ACC_Aclm_OR_Pzplt_indv,    ACC_Acstc_OR_Pzplt_indv,   ACC_all_sensors_OR_cmbd_indv, ...
    ACC_Aclm_AND_Acstc_indv,   ACC_Aclm_AND_Pzplt_indv,   ACC_Acstc_AND_Pzplt_indv,  ACC_all_sensors_AND_cmbd_indv, ...
    ACC_atleast_1_sensor_indv, ACC_atleast_2_sensor_indv, ACC_atleast_3_sensor_indv, ACC_atleast_4_sensor_indv,     ACC_atleast_5_sensor_indv, ACC_atleast_6_sensor_indv, ...
    ACC_atleast_1_type_indv,   ACC_atleast_2_type_indv,   ACC_atleast_3_type_indv];
TPD_indv_total  = [TPD_all_indv_sensors_indv, TPD_Aclm_indv, TPD_Acstc_indv, TPD_Pzplt_indv,...
    TPD_Aclm_OR_Acstc_indv,    TPD_Aclm_OR_Pzplt_indv,    TPD_Acstc_OR_Pzplt_indv,   TPD_all_sensors_OR_cmbd_indv, ...
    TPD_Aclm_AND_Acstc_indv,   TPD_Aclm_AND_Pzplt_indv,   TPD_Acstc_AND_Pzplt_indv,  TPD_all_sensors_AND_cmbd_indv, ...
    TPD_atleast_1_sensor_indv, TPD_atleast_2_sensor_indv, TPD_atleast_3_sensor_indv, TPD_atleast_4_sensor_indv,     TPD_atleast_5_sensor_indv, TPD_atleast_6_sensor_indv, ...
    TPD_atleast_1_type_indv,   TPD_atleast_2_type_indv,   TPD_atleast_3_type_indv];
FPD_indv_total  = [FPD_all_indv_sensors_indv, FPD_Aclm_indv, FPD_Acstc_indv, FPD_Pzplt_indv,...
    FPD_Aclm_OR_Acstc_indv,    FPD_Aclm_OR_Pzplt_indv,    FPD_Acstc_OR_Pzplt_indv,   FPD_all_sensors_OR_cmbd_indv, ...
    FPD_Aclm_AND_Acstc_indv,   FPD_Aclm_AND_Pzplt_indv,   FPD_Acstc_AND_Pzplt_indv,  FPD_all_sensors_AND_cmbd_indv, ...
    FPD_atleast_1_sensor_indv, FPD_atleast_2_sensor_indv, FPD_atleast_3_sensor_indv, FPD_atleast_4_sensor_indv,     FPD_atleast_5_sensor_indv, FPD_atleast_6_sensor_indv, ...
    FPD_atleast_1_type_indv,   FPD_atleast_2_type_indv,   FPD_atleast_3_type_indv];
TND_indv_total  = [TND_all_indv_sensors_indv, TND_Aclm_indv, TND_Acstc_indv, TND_Pzplt_indv,...
    TND_Aclm_OR_Acstc_indv,    TND_Aclm_OR_Pzplt_indv,    TND_Acstc_OR_Pzplt_indv,   TND_all_sensors_OR_cmbd_indv, ...
    TND_Aclm_AND_Acstc_indv,   TND_Aclm_AND_Pzplt_indv,   TND_Acstc_AND_Pzplt_indv,  TND_all_sensors_AND_cmbd_indv, ...
    TND_atleast_1_sensor_indv, TND_atleast_2_sensor_indv, TND_atleast_3_sensor_indv, TND_atleast_4_sensor_indv,     TND_atleast_5_sensor_indv, TND_atleast_6_sensor_indv, ...
    TND_atleast_1_type_indv,   TND_atleast_2_type_indv,   TND_atleast_3_type_indv];
FND_indv_total  = [FND_all_indv_sensors_indv, FND_Aclm_indv, FND_Acstc_indv, FND_Pzplt_indv,...
    FND_Aclm_OR_Acstc_indv,    FND_Aclm_OR_Pzplt_indv,    FND_Acstc_OR_Pzplt_indv,   FND_all_sensors_OR_cmbd_indv, ...
    FND_Aclm_AND_Acstc_indv,   FND_Aclm_AND_Pzplt_indv,   FND_Acstc_AND_Pzplt_indv,  FND_all_sensors_AND_cmbd_indv, ...
    FND_atleast_1_sensor_indv, FND_atleast_2_sensor_indv, FND_atleast_3_sensor_indv, FND_atleast_4_sensor_indv,     FND_atleast_5_sensor_indv, FND_atleast_6_sensor_indv, ...
    FND_atleast_1_type_indv,   FND_atleast_2_type_indv,   FND_atleast_3_type_indv];

FS_overall_total  = [FS_all_indv_sensors_overall, FS_Aclm_overall, FS_Acstc_overall, FS_Pzplt_overall,...
    FS_Aclm_OR_Acstc_overall,    FS_Aclm_OR_Pzplt_overall,    FS_Acstc_OR_Pzplt_overall,   FS_all_sensors_OR_cmbd_overall, ...
    FS_Aclm_AND_Acstc_overall,   FS_Aclm_AND_Pzplt_overall,   FS_Acstc_AND_Pzplt_overall,  FS_all_sensors_AND_cmbd_overall, ...
    FS_atleast_1_sensor_overall, FS_atleast_2_sensor_overall, FS_atleast_3_sensor_overall, FS_atleast_4_sensor_overall,     FS_atleast_5_sensor_overall, FS_atleast_6_sensor_overall, ...
    FS_atleast_1_type_overall,   FS_atleast_2_type_overall,   FS_atleast_3_type_overall];
SEN_overall_total  = [SEN_all_indv_sensors_overall, SEN_Aclm_overall, SEN_Acstc_overall, SEN_Pzplt_overall,...
    SEN_Aclm_OR_Acstc_overall,    SEN_Aclm_OR_Pzplt_overall,    SEN_Acstc_OR_Pzplt_overall,   SEN_all_sensors_OR_cmbd_overall, ...
    SEN_Aclm_AND_Acstc_overall,   SEN_Aclm_AND_Pzplt_overall,   SEN_Acstc_AND_Pzplt_overall,  SEN_all_sensors_AND_cmbd_overall, ...
    SEN_atleast_1_sensor_overall, SEN_atleast_2_sensor_overall, SEN_atleast_3_sensor_overall, SEN_atleast_4_sensor_overall,     SEN_atleast_5_sensor_overall, SEN_atleast_6_sensor_overall, ...
    SEN_atleast_1_type_overall,   SEN_atleast_2_type_overall,   SEN_atleast_3_type_overall];
PPV_overall_total  = [PPV_all_indv_sensors_overall, PPV_Aclm_overall, PPV_Acstc_overall, PPV_Pzplt_overall,...
    PPV_Aclm_OR_Acstc_overall,    PPV_Aclm_OR_Pzplt_overall,    PPV_Acstc_OR_Pzplt_overall,   PPV_all_sensors_OR_cmbd_overall, ...
    PPV_Aclm_AND_Acstc_overall,   PPV_Aclm_AND_Pzplt_overall,   PPV_Acstc_AND_Pzplt_overall,  PPV_all_sensors_AND_cmbd_overall, ...
    PPV_atleast_1_sensor_overall, PPV_atleast_2_sensor_overall, PPV_atleast_3_sensor_overall, PPV_atleast_4_sensor_overall,     PPV_atleast_5_sensor_overall, PPV_atleast_6_sensor_overall, ...
    PPV_atleast_1_type_overall,   PPV_atleast_2_type_overall,   PPV_atleast_3_type_overall];
SPE_overall_total  = [SPE_all_indv_sensors_overall, SPE_Aclm_overall, SPE_Acstc_overall, SPE_Pzplt_overall,...
    SPE_Aclm_OR_Acstc_overall,    SPE_Aclm_OR_Pzplt_overall,    SPE_Acstc_OR_Pzplt_overall,   SPE_all_sensors_OR_cmbd_overall, ...
    SPE_Aclm_AND_Acstc_overall,   SPE_Aclm_AND_Pzplt_overall,   SPE_Acstc_AND_Pzplt_overall,  SPE_all_sensors_AND_cmbd_overall, ...
    SPE_atleast_1_sensor_overall, SPE_atleast_2_sensor_overall, SPE_atleast_3_sensor_overall, SPE_atleast_4_sensor_overall,     SPE_atleast_5_sensor_overall, SPE_atleast_6_sensor_overall, ...
    SPE_atleast_1_type_overall,   SPE_atleast_2_type_overall,   SPE_atleast_3_type_overall];
ACC_overall_total  = [ACC_all_indv_sensors_overall, ACC_Aclm_overall, ACC_Acstc_overall, ACC_Pzplt_overall,...
    ACC_Aclm_OR_Acstc_overall,    ACC_Aclm_OR_Pzplt_overall,    ACC_Acstc_OR_Pzplt_overall,   ACC_all_sensors_OR_cmbd_overall, ...
    ACC_Aclm_AND_Acstc_overall,   ACC_Aclm_AND_Pzplt_overall,   ACC_Acstc_AND_Pzplt_overall,  ACC_all_sensors_AND_cmbd_overall, ...
    ACC_atleast_1_sensor_overall, ACC_atleast_2_sensor_overall, ACC_atleast_3_sensor_overall, ACC_atleast_4_sensor_overall,     ACC_atleast_5_sensor_overall, ACC_atleast_6_sensor_overall, ...
    ACC_atleast_1_type_overall,   ACC_atleast_2_type_overall,   ACC_atleast_3_type_overall];
TPD_overall_total  = [TPD_all_indv_sensors_overall, TPD_Aclm_overall, TPD_Acstc_overall, TPD_Pzplt_overall,...
    TPD_Aclm_OR_Acstc_overall,    TPD_Aclm_OR_Pzplt_overall,    TPD_Acstc_OR_Pzplt_overall,   TPD_all_sensors_OR_cmbd_overall, ...
    TPD_Aclm_AND_Acstc_overall,   TPD_Aclm_AND_Pzplt_overall,   TPD_Acstc_AND_Pzplt_overall,  TPD_all_sensors_AND_cmbd_overall, ...
    TPD_atleast_1_sensor_overall, TPD_atleast_2_sensor_overall, TPD_atleast_3_sensor_overall, TPD_atleast_4_sensor_overall,     TPD_atleast_5_sensor_overall, TPD_atleast_6_sensor_overall, ...
    TPD_atleast_1_type_overall,   TPD_atleast_2_type_overall,   TPD_atleast_3_type_overall];
FPD_overall_total  = [FPD_all_indv_sensors_overall, FPD_Aclm_overall, FPD_Acstc_overall, FPD_Pzplt_overall,...
    FPD_Aclm_OR_Acstc_overall,    FPD_Aclm_OR_Pzplt_overall,    FPD_Acstc_OR_Pzplt_overall,   FPD_all_sensors_OR_cmbd_overall, ...
    FPD_Aclm_AND_Acstc_overall,   FPD_Aclm_AND_Pzplt_overall,   FPD_Acstc_AND_Pzplt_overall,  FPD_all_sensors_AND_cmbd_overall, ...
    FPD_atleast_1_sensor_overall, FPD_atleast_2_sensor_overall, FPD_atleast_3_sensor_overall, FPD_atleast_4_sensor_overall,     FPD_atleast_5_sensor_overall, FPD_atleast_6_sensor_overall, ...
    FPD_atleast_1_type_overall,   FPD_atleast_2_type_overall,   FPD_atleast_3_type_overall];
TND_overall_total  = [TND_all_indv_sensors_overall, TND_Aclm_overall, TND_Acstc_overall, TND_Pzplt_overall,...
    TND_Aclm_OR_Acstc_overall,    TND_Aclm_OR_Pzplt_overall,    TND_Acstc_OR_Pzplt_overall,   TND_all_sensors_OR_cmbd_overall, ...
    TND_Aclm_AND_Acstc_overall,   TND_Aclm_AND_Pzplt_overall,   TND_Acstc_AND_Pzplt_overall,  TND_all_sensors_AND_cmbd_overall, ...
    TND_atleast_1_sensor_overall, TND_atleast_2_sensor_overall, TND_atleast_3_sensor_overall, TND_atleast_4_sensor_overall,     TND_atleast_5_sensor_overall, TND_atleast_6_sensor_overall, ...
    TND_atleast_1_type_overall,   TND_atleast_2_type_overall,   TND_atleast_3_type_overall];
FND_overall_total  = [FND_all_indv_sensors_overall, FND_Aclm_overall, FND_Acstc_overall, FND_Pzplt_overall,...
    FND_Aclm_OR_Acstc_overall,    FND_Aclm_OR_Pzplt_overall,    FND_Acstc_OR_Pzplt_overall,   FND_all_sensors_OR_cmbd_overall, ...
    FND_Aclm_AND_Acstc_overall,   FND_Aclm_AND_Pzplt_overall,   FND_Acstc_AND_Pzplt_overall,  FND_all_sensors_AND_cmbd_overall, ...
    FND_atleast_1_sensor_overall, FND_atleast_2_sensor_overall, FND_atleast_3_sensor_overall, FND_atleast_4_sensor_overall,     FND_atleast_5_sensor_overall, FND_atleast_6_sensor_overall, ...
    FND_atleast_1_type_overall,   FND_atleast_2_type_overall,   FND_atleast_3_type_overall];

% Initialization of the table
T_sensor_combinations_overall = table('Size',[length(sensor_combinations) 13], 'VariableTypes', {'string', 'double', 'double', 'string', 'double', ...
    'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'}) ; 
T_sensor_combinations_indv = table('Size',[length(sensor_combinations)*n_data_files 14], 'VariableTypes', {'string', 'string', 'double', 'double', ...
    'string', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'}) ; 


% Loop for storing the data
for i = 1 : length(sensor_combinations)

    T_sensor_combinations_overall{i,:} = [sensor_combinations(i), sum(duration_data_files), mean(Force_mean), FM_min_SN_used, FS_overall_total(i), ...
        SEN_overall_total(i), PPV_overall_total(i), SPE_overall_total(i), ACC_overall_total(i), TPD_overall_total(i), FPD_overall_total(i), ...
        TND_overall_total(i), FND_overall_total(i)];

    for j = 1 : n_data_files
        if ischar(data_file_names) % If there is only a single data file
            DFN = data_file_names;
        else
            DFN = cellstr(data_file_names{j}); % Data file names are converted from cell elements to strings
        end

        T_sensor_combinations_indv{(i-1)*n_data_files + j, :} = [sensor_combinations(i), DFN, duration_data_files(j), Force_mean(j), ...
            FM_min_SN_used, FS_indv_total{i}(j), SEN_indv_total{i}(j), PPV_indv_total{i}(j), SPE_indv_total{i}(j), ACC_indv_total{i}(j), ...
            TPD_indv_total{i}(j), FPD_indv_total{i}(j), TND_indv_total{i}(j), FND_indv_total{i}(j)];
    end

end

T_sensor_combinations_overall.Properties.VariableNames = {'Sensor Name', 'Duration (min)', 'Mean Force (au)', 'Threshold multiplier', 'F Score', ...
    'Sensitivity', 'PPV', 'Specificity', 'Accuracy', 'TPD', 'FPD', 'TND', 'FND'}; % Assigns column names to the table
T_sensor_combinations_indv.Properties.VariableNames = {'Sensor Name', 'Data file', 'Duration (min)', 'Mean Force (au)', 'Threshold multiplier', 'F Score', ...
    'Sensitivity', 'PPV', 'Specificity', 'Accuracy', 'TPD', 'FPD', 'TND', 'FND'}; % Assigns column names to the table
%

% --------------------- Displaying the stored data ---------------------- %
% disp('Performance of overall data sets: ');
T_sensor_combinations_final=T_sensor_combinations_overall([7 8 9 13 17 25], :);
disp(T_sensor_combinations_final)
%

%% EXTRACTION OF DETECTIONS FOR T-F ANALYSIS ==============================
% TPDs and FPDs that are common to all the sensors are extracted for T-F
% analysis

% Starting notification
disp('T-F analysis is going on... ')

% Parameters and variables for segmentation and detection matching
ext_backward = 5.0; % Backward extension length in second
ext_forward = 2.0; % Forward extension length in second
FM_dilation_time = 3.0;% Dilation size in seconds
n_FM_sensors = 6; % number of FM sensors
FM_min_SN = [30, 30, 60, 60, 50, 50];

IMU_map = cell(1,n_data_files);   
M_sntn_map = cell(1,n_data_files);  
threshold = zeros(n_data_files, n_FM_sensors);

% Parameters and variables for PSD estimation
w = hann(Fs_sensor); % Hann window of length 1s
n_overlap = Fs_sensor/2; % 50% overlap between the segments
p_overlap = 0.5; % 50% overlap between the segments
nfft = Fs_sensor; % number of fourier points will be same as sampling freq

Pxx_TPD_avg = zeros(Fs_sensor/2 + 1, n_FM_sensors);
Pxx_detrended_TPD_avg = zeros(Fs_sensor/2 + 1, n_FM_sensors);
Pxx_FPD_avg = zeros(Fs_sensor/2 + 1, n_FM_sensors);
Pxx_detrended_FPD_avg = zeros(Fs_sensor/2 + 1, n_FM_sensors);
n_TPD_cycle = 0;
n_FPD_cycle = 0;

TPD_extracted = cell(1, n_data_files);
FPD_extracted = cell(1, n_data_files);
 
for i = 1 : n_data_files

    % Starting notification
    fprintf('Current data file: %.0f/%.0f\n', i,n_data_files)

    % --------- Segmentaiton of IMU data and creation of IMU_map ---------%
    % get_IMU_map() function is used here. It segments and dilates the IMU 
    % data and returns the resultant data as IMU map. Settings for the
    % segmentation and dilation are given inside the function.
    %   Input variables:  IMU_data- A column vector of IMU data 
    %                     data_file_names- a char variable with data file name
    %   Output variables: IMU_map- A column vector with segmented IMU data 

    if n_data_files == 1
        IMU_map{i} = get_IMU_map(IMU_data_fltd{i}, data_file_names, Fs_sensor);
    else
        IMU_map{i} = get_IMU_map(IMU_data_fltd{i}, data_file_names{i}, Fs_sensor);
    end

    % ----------------------- Creation of M_sensation_map ----------------%
    % get_sensation_map() function is used here. This function gives the
    % sensation map by dilating every detection to past and future. It also 
    % revomes the windows that overlaps with IMU_map. Settings for the 
    % extension are given inside the function.
    %   Input variables-  sensation_data- a column vector with the sensation data
    %                     IMU_map- a column vector
    %                     ext_backward, ext_forward- scalar values
    %                     Fs_sensor,Fs_sensation- scalar values    %                     
    %   Output variables: M_sntn_map- a column vector with all the sensation windows joined together 

    M_sntn_map{i} = get_sensation_map(sensation_data_SD1_trimd{i}, IMU_map{i}, ext_backward, ext_forward, Fs_sensor, Fs_sensation);

    % ~~~~~~~~~~~~~~~~~~~~~~ Segmentation of FM data ~~~~~~~~~~~~~~~~~~~~~%
    % get_segmented_data() function will be used here, which will
    % threshold the data, remove body movement, and dilate the data.
    % Setting for the threshold and dilation are given in the function.
    %   Input variables:  sensor_data- a cell variable;
    %                     min_SN- a vector/ a scalar;
    %                     IMU_map- a column vector
    %                     Fs_sensor, FM_dilation_time- a scalar
    %   Output variables: sensor_data_sgmntd- a cell variable of same size as the input variable sensor_data_fltd;
    %                     h- a vector

    sensor_data_fltd = {Aclm_data1_fltd{i}, Aclm_data2_fltd{i}, Acstc_data1_fltd{i}, Acstc_data2_fltd{i}, ...
        Pzplt_data1_fltd{i}, Pzplt_data2_fltd{i}};
    [sensor_data_sgmntd, threshold(i,:)] = get_segmented_data(sensor_data_fltd, FM_min_SN, IMU_map{i}, FM_dilation_time, Fs_sensor);

    % --------------- Detections common to all the sensors ----------------
    sensor_data_sgmntd_cmbd_all_sensors = {double(sensor_data_sgmntd{1} | sensor_data_sgmntd{2} | sensor_data_sgmntd{3} | ...
        sensor_data_sgmntd{4} | sensor_data_sgmntd{5} | sensor_data_sgmntd{6})};
    sensor_data_sgmntd_cmbd_all_sensors_labeled = bwlabel(sensor_data_sgmntd_cmbd_all_sensors{1});
    n_label = length(unique(sensor_data_sgmntd_cmbd_all_sensors_labeled)) - 1; % Number of labels in the sensor_data_cmbd_all_sensors_labeled

    % Initialization of variables
    sensor_data_sgmntd_atleast_6_sensor{1} = zeros(length(sensor_data_sgmntd{1}),1);

    if (n_label) % When there is a detection by the sensor system

        for k = 1 : n_label
            L_min = find(sensor_data_sgmntd_cmbd_all_sensors_labeled == k, 1 ); % Sample no. corresponding to start of the label
            L_max = find(sensor_data_sgmntd_cmbd_all_sensors_labeled == k, 1, 'last' ); % Sample no. corresponding to end of the label

            indv_detection_map = zeros(length(sensor_data_sgmntd{1}),1); % Need to be initialized before every detection matching
            indv_detection_map(L_min:L_max) = 1; % mapping individual sensation data

            % For detection by at least n sensors
            tmp_var = 0; % Variable to hold number of common sensors for each detection
            for j = 1:n_FM_sensors
                if (sum(indv_detection_map.*sensor_data_sgmntd{j})) % Non-zero value indicates intersection
                    tmp_var = tmp_var + 1;
                end
            end

            switch tmp_var
                case 6
                    sensor_data_sgmntd_atleast_6_sensor{1}(L_min:L_max) = 1;
            end
        end
    end

    % ------------------- Extraction of common TPDs and FPDs---------------
    sensor_data_sgmntd_atleast_6_sensor_labeled = bwlabel(sensor_data_sgmntd_atleast_6_sensor{1});
    n_detection = length(unique(sensor_data_sgmntd_atleast_6_sensor_labeled)) - 1; % Number of total detections 
    n_candidate_TPD = length(unique(sensor_data_sgmntd_atleast_6_sensor_labeled.*M_sntn_map{i})) - 1; % Number of detections that intersects with maternal sensation
    n_candidate_FPD =  n_detection - n_candidate_TPD;

    current_file_TPD_extraction_cell = cell(1, n_candidate_TPD); % Each cell will contain 1 TPD for all 6 sensors
    current_file_FPD_extraction_cell = cell(1, n_candidate_FPD);

    if n_detection
        k_TPD = 0;
        k_FPD = 0;
        for k = 1:n_detection
            L_min = find(sensor_data_sgmntd_atleast_6_sensor_labeled == k, 1 ); % Sample no. corresponding to the start of the label
            L_max = find(sensor_data_sgmntd_atleast_6_sensor_labeled == k, 1, 'last' ); % Sample no. corresponding to the end of the label
            indv_window = zeros(length(M_sntn_map{i}),1); % This window represents the current detection, which will be compared with M_sntn_map to find if it is TPD or FPD
            indv_window(L_min:L_max) = 1;
            X = sum(indv_window.*M_sntn_map{i}); % Checks the overlap with the maternal sensation

            if X
                k_TPD = k_TPD + 1;
                current_TPD_extraction = zeros(L_max-L_min+1, n_FM_sensors);
                for j = 1:n_FM_sensors
                    current_TPD_extraction(:,j) = sensor_data_fltd{j}(L_min:L_max); % Each 

                    [Pxx_TPD, f_TPD] = pwelch(current_TPD_extraction(:,j), w, n_overlap, nfft, Fs_sensor);
                    [Pxx_detrended_TPD, f_detrended_TPD] = pwelch_new(current_TPD_extraction(:,j), w, p_overlap, nfft, Fs_sensor, 'half', 'plot', 'mean');
                    % pwelch_new() is a user defined function that allows detrending.
                    %   Input conditions: half - one sided PSD; plot- normal plotting; 
                    %                     mean-  remove the mean value of each segment from each segment of the data.

                    Pxx_TPD_avg(:,j) = Pxx_TPD_avg(:,j) + Pxx_TPD; % Summing the PSD for each TPD
                    Pxx_detrended_TPD_avg(:,j) = Pxx_detrended_TPD_avg(:,j) + Pxx_detrended_TPD;
                end
                current_file_TPD_extraction_cell {k_TPD} = current_TPD_extraction;
                n_TPD_cycle = n_TPD_cycle + 1;
            else
                k_FPD = k_FPD + 1;
                current_FPD_extraction = zeros(L_max-L_min+1, n_FM_sensors);
                for j = 1:n_FM_sensors
                    current_FPD_extraction(:,j) = sensor_data_fltd{j}(L_min:L_max);

                    [Pxx_FPD, f_FPD] = pwelch(current_FPD_extraction(:,j), w, n_overlap, nfft, Fs_sensor);
                    [Pxx_detrended_FPD, f_detrended_FPD] = pwelch_new(current_FPD_extraction(:,j), w, p_overlap, nfft, Fs_sensor, 'half', 'plot', 'mean');

                    Pxx_FPD_avg(:,j) = Pxx_FPD_avg(:,j) + Pxx_FPD;
                    Pxx_detrended_FPD_avg(:,j) = Pxx_detrended_FPD_avg(:,j) + Pxx_detrended_FPD;
                end
                current_file_FPD_extraction_cell {k_FPD} = current_FPD_extraction;
                n_FPD_cycle = n_FPD_cycle + 1;
            end
        end
    end

    TPD_extracted{i} = current_file_TPD_extraction_cell; % Each cell will contain wll the extracted TPD in a particular data set
    FPD_extracted{i} = current_file_FPD_extraction_cell;
end

% Averaging the PSD summation
Pxx_TPD_avg = Pxx_TPD_avg/n_TPD_cycle; 
Pxx_detrended_TPD_avg =  Pxx_detrended_TPD_avg/n_TPD_cycle;
Pxx_FPD_avg = Pxx_FPD_avg/n_FPD_cycle;
Pxx_detrended_FPD_avg =  Pxx_detrended_FPD_avg/n_FPD_cycle;

%% SPECTROGRAM ============================================================

% Parameters and variables for PSD estimation
s_size = Fs_sensor/2; % Size of each segment in the STFT
w = hann(s_size); % Hann window of length 0.5 S
n_overlap = floor(s_size/1.25); % 80% overlap between the segments
nfft = s_size*2; % number of fourier points will be twice sampling freq. This pads 0 at the end to increase the frequency resolution
% For the above setting, the obtained time resolution = 100 ms, and frequency resolution = 1 Hz 

% Selection of data file and TPD
data_file_no = 5;
TPD_no = 7;
t_start = 2; % Starting time in s
t_end = 7; % Ending time in s

% Good candidates: (File_no, TPD_no): (3,3), (4,2), (4,3) , (5,2), (5,4), (5,5)
% excellent,  (5,6, excellent^2), (5,7, excellent^3), (5,11, too good),
% (6,6, too good), (7,2, super)
% Finally used: (5,7)

% Plot settings
legend_all={'Accelerometer','Acoustic sensor','Piezoelectric diaphragm'};
L_width = 3; % Width of the box
p_L_width = 4; % Width of the lines in the plot
F_size = 28; % Font size
F_name = 'Times New Roman';

tiledlayout(3,3,'Padding','tight','TileSpacing','tight'); % 3x3 tile in most compact fitting

% Plot the data in the first row
for j = 1:3
    % Get the X and Y-axis data
    data = TPD_extracted{1, data_file_no}{1, TPD_no}(:, 2*j); % TO get the data for the right sensor only
    % data = data(t_start*Fs_sensor:t_end*Fs_sensor+1,1); % uncomment if you want to plot within a range defined by t_start:t_end
    time_vec = (0:length(data)-1)/Fs_sensor;

    % Plot the data
    nexttile
    plot(time_vec, data, 'LineWidth', p_L_width)    

    % Set the axis properties
    ax = gca; % To get the coordinate axis
    if j == 1
        ylabel('Amplitude (a.u.)'); % Ylabel is only applied to the leftmost axis
        ax.YAxis.Exponent = -2;
        ylim([-0.02 0.01])
    end
    if j == 2
        ax.YAxis.Exponent = -2;
        ylim([-0.06 0.06])
        yticks(-0.06:0.03:0.06)
    end
    if j == 3
        ax.YAxis.Exponent = -1;
        ylim([-0.3 0.3])
        yticks(-0.3:0.1:0.3)
    end

    xlabel('Time(s)');
    xlim([0 8])
    xticks(0:2:10);
    set(gca, 'FontName', F_name, 'FontSize', F_size, 'linewidth', L_width) % Sets the font type and size of the labels
end

% Plot the spectrogram in the seond row
for j=1:3   
    % Get the X and Y-axis data
    data = TPD_extracted{1,data_file_no}{1,TPD_no}(:,2*j);
    % data = data(t_start*Fs_sensor:t_end*Fs_sensor+1,1); % Uncomment if you want to plot within a range defined by t_start:t_end
    [~,f,t,P] = spectrogram(data, w, n_overlap, nfft, Fs_sensor, 'psd', 'yaxis'); % P is the PSD 

    % Plot the data
    nexttile
    % spectrogram(data,w_spec,n_overlap_spec,nfft_spec,Fs_sensor,'yaxis') % Plot using spectrogram
    imagesc(t, f, (P+eps)) % Plot using imagesc. 
    %      Add eps like pspectrogram does: eps is the distance between 1.0 and the next largest double precision number (eps = 2.2204e-16)
    axis xy % Corrects the Y-axis order: By default imagesc() puts Y-axis to high-to-low order
    h = colorbar; % To display the colorbar on the left    

    % Set the axis properties
    if j == 3
        h.Label.String = 'PSD (a.u.^2/Hz)'; % Labels the colorbar
        h.Ticks = 3*10^-3:3*10^-3:12*10^-3  ; % To manage colorbar ticks
    end
    colormap(jet)
    ylim([0,30]);
    yticks(0:10:30);
    if j==1
        ylabel('Frequency (Hz)')
    end
    xlim([-inf 8])
    xlabel('Time (s)'); 
    set(gca, 'FontName', F_name, 'FontSize', F_size, 'linewidth', L_width) % Sets the font type & size, and the line width
end

% Plot the PSD in the third row
for j=1:3
    % Plot the data
    nexttile
    plot(f_detrended_FPD, Pxx_detrended_FPD_avg(:,2*j),'LineWidth', p_L_width)

    % Set the axis properties
    if j == 1
        ylabel('PSD(a.u.^2/Hz)')
        ylim([0 8*10^-7])
        yticks(0:2*10^-7:8*10^-7)
    end
    if j == 2
        ylim([0 3*10^-6])
        yticks(0:1*10^-6:3*10^-6)
    end
    if j == 3
        ylim([0 8*10^-5])
        yticks(0:2*10^-5:8*10^-5)
    end

    xlim([0 30])
    xticks(0:10:30)
    xlabel('Frequency (Hz)')
    set(gca, 'FontName', F_name, 'FontSize', F_size, 'linewidth', L_width) % Sets the font type & size, and the line width
end
%
