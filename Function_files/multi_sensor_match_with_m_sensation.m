function [TPD, FPD, TND, FND, Detected_M_sntn_Map, Detected_TPD_Map_indv_sensor, sensor_data_sgmntd_removed, sensor_data_fltd_removed] = multi_sensor_match_with_m_sensation(sensor_data_fltd, ...
    sensor_data_sgmntd, sensation_data, IMU_map, M_sntn_Map, ext_backward, ext_forward, FM_dilation_time, Fs_sensor, Fs_sensation, rqrd_cmn)

%MATCH_WITH_M_SENSATION Summary of this function goes here
%   Input variables:  sensor_data_fltd, sensor_data_sgmntd- Cell matrices with data for each sensor in each cell. Each cell are column vectors
%                     sensation_data, IMU_map, M_sntn_Map- Cell matrices with data for each file in each cell. Each cell is a colunm vector.
%                     ext_backward, ext_forward, FM_dilation_time- scalar values
%                     Fs_sensor, Fs_sensation, rqrd_cmn- scalar values
%
%   Output variables: TPD, FPD, TND, FND- scalar values 
%                     Detected_M_sntn_Map, Detected_TPD_Map_indv_sensor-
%                     sensor_data_sgmntd_removed, sensor_data_fltd_removed-                      
% 

% Calculation of fixed values
n_sensors = length(sensor_data_sgmntd); % sensor_data_sgmntd is cell variable with data for each sensor in each cell
matching_window_size = ext_backward + ext_forward; % window size is equal to the window size used to create the maternal sensation map
minm_overlap_time = FM_dilation_time/2; % Minimum overlap in second
DLB = round(ext_backward*Fs_sensor); % backward extension length in sample number
DLF = round(ext_forward*Fs_sensor); % forward extension length in sample number

% Labeling sensation data and determining number of maternal sensation detection
sensation_data_labeled = bwlabel(sensation_data);
n_Maternal_detected_movement = length(unique(sensation_data_labeled)) - 1; % 1 is deducted to remove the initial value of 0

% Variable decleration
TPD_all = zeros(n_Maternal_detected_movement, n_sensors); % True positive detection
FND_all = zeros(n_Maternal_detected_movement, n_sensors); % False negative detection

% Loop for matching individual maternal sensation detection
for j = 1 : n_sensors

    % ------------------ Determination of TPD and FND ----------------%
    if (n_Maternal_detected_movement) % When there is a detection by the mother

        for k = 1 : n_Maternal_detected_movement

            L_min = find(sensation_data_labeled == k, 1 ); % Sample no. corresponding to start of the label
            L_max = find(sensation_data_labeled == k, 1, 'last' ); % Sample no. corresponding to end of the label

            L1 = L_min*round(Fs_sensor/Fs_sensation) - DLB; % sample no. for the starting point of this sensation in the map
            L1 = max(L1,1); % Just a check so that L1 remains higher than 1st data sample

            L2 = L_max*round(Fs_sensor/Fs_sensation) + DLF; % sample no. for the ending point of this sensation in the map
            L2 = min(L2,length(sensor_data_sgmntd{j})); % Just a check so that L2 remains lower than the last data sample

            indv_sensation_map = zeros(length(sensor_data_sgmntd{j}),1); % Need to be initialized before every detection matching
            indv_sensation_map(L1:L2) = 1; % mapping individual sensation data

            X = sum(indv_sensation_map.*IMU_map); % this is non-zero if there is a coincidence with maternal body movement

            if (~X)  % true when there is no coincidence with body movement
                % TPD and FND calculation for indivudial sensors
                Y = sum(sensor_data_sgmntd{j}.*indv_sensation_map); % Non-zero value gives the matching
                if (Y) % true if there is a coincidence
                    TPD_all(k,j) = 1; % TPD assigned
                else
                    FND_all(k,j) = 1; % FND assigned
                end
            end

        end

    end
end

% Calculation of TPD and FND for combined sensors
TPD = sum(sum(TPD_all,2) >= rqrd_cmn); % summation along the rows
FND = sum(sum(FND_all,2) > (n_sensors - rqrd_cmn));
%

% --------------------- Determination of TND and FPD  --------------------%
% Creation of a map for the detected region
detected_m_sntn = sum(TPD_all,2) >= rqrd_cmn; % Gives the location of detected M_sntns as 1
index_detected_m_sntn = find(detected_m_sntn); % Gives the index of detected M_sntns in terms of sensation number in sensation_data_labeled
n_detected_m_sntn = length(index_detected_m_sntn); % Number of detected M_sntn; should be same as TPD

Detected_M_sntn_Map = zeros(length(M_sntn_Map),1);  % Initialization of the map that will contain only detected M_sntn

if (n_detected_m_sntn) % When there is a TPD by the sensor

    for k = 1 : n_detected_m_sntn

        L_min = find(sensation_data_labeled == index_detected_m_sntn(k), 1 ); % Sample no. corresponding to start of the label
        L_max = find(sensation_data_labeled == index_detected_m_sntn(k), 1, 'last' ); % Sample no. corresponding to end of the label

        L1 = L_min*round(Fs_sensor/Fs_sensation) - DLB; % sample no. for the starting point of this sensation in the map
        L1 = max(L1,1); % Just a check so that L1 remains higher than 1st data sample

        L2 = L_max*round(Fs_sensor/Fs_sensation) + DLF; % sample no. for the ending point of this sensation in the map
        L2 = min(L2,length(M_sntn_Map)); % Just a check so that L2 remains lower than the last data sample

        Detected_M_sntn_Map(L1:L2) = 1; % mapping of detected sensation data

    end

end
%

% Creation of a map that will hold the TPDs from all the sensors
Detected_TPD_Map_all_sensor = zeros(length(M_sntn_Map),1); % This map will hold the TPD region from all the sensors
Detected_TPD_Map_indv_sensor = zeros(length(M_sntn_Map),n_sensors); % This map will hold the TPD regions for individual sensors in each column

for j = 1 : n_sensors

    sensor_data_sgmntd_labeled = bwlabel(sensor_data_sgmntd{j});
    curnt_matched_vector = sensor_data_sgmntd_labeled.*Detected_M_sntn_Map; % Non-zero elements gives the matching
    curnt_matched_label = unique(curnt_matched_vector); % Gives the label of the matched sensor data segments

    if (length(curnt_matched_label)>1)
        curnt_matched_label = curnt_matched_label(2:end); % Removes the first element, which is 0

        for m = 1 : length(curnt_matched_label)
            Detected_TPD_Map_all_sensor(sensor_data_sgmntd_labeled == curnt_matched_label(m)) = 1; % Assigns 1 to the TPD segments of the segmented signal
            Detected_TPD_Map_indv_sensor(sensor_data_sgmntd_labeled == curnt_matched_label(m),j) = 1; % Holds the detected TPD for each sensors in each row
        end

    end

end
%

% Loop for detecting FPD and TND
% sensor_data_sgmntd_removed = cell(1, n_sensors); % Varaible to store segmented data after removing the TPD and M_sensation regions
for j = 1 : n_sensors

    % Removal of the TPD and FND parts from the individual sensor data
    Arb_val = 4; % An arbitrary value
    sensor_data_sgmntd{j}(Detected_TPD_Map_all_sensor == 1) = Arb_val; % Assigns an arbitrary value to the TPD areas
    sensor_data_sgmntd{j}(M_sntn_Map == 1) = Arb_val; % Assigns an arbitrary value to the area under the M_sntn_Map
    sensor_data_sgmntd_removed (:,j) = sensor_data_sgmntd{j}(sensor_data_sgmntd{j} ~= Arb_val); % Removes all the elements with value = Arb_val from the segmented data
    sensor_data_fltd_removed (:,j) = sensor_data_fltd{j}(sensor_data_sgmntd{j} ~= Arb_val); 
    %

    % Calculation of TNDs and FPDs for individual sensors
    L_removed = length(sensor_data_sgmntd_removed);
    index_window_start = 1;
    index_window_end = min(index_window_start + Fs_sensor*matching_window_size, L_removed);
    i = 1; % index of epoches

    while (index_window_start < L_removed)
        indv_window = sensor_data_sgmntd_removed(index_window_start: index_window_end, j);
        index_non_zero = find(indv_window);

        if length(index_non_zero) >= (minm_overlap_time*Fs_sensor)
            FPD_all(i,j) = 1;
        else
            TND_all(i,j) = 1;
        end

        i = i+1;
        index_window_start = index_window_end + 1;
        index_window_end = min(index_window_start + Fs_sensor*matching_window_size, L_removed);
    end

end

% Calculation of FPD and TND for the combined sensors
FPD = sum(sum(FPD_all,2) >= rqrd_cmn); % summation along the rows
TND = sum(sum(TND_all,2) > (n_sensors - rqrd_cmn));

end

