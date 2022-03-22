function [TPD, FPD, TND, FND] = match_with_m_sensation(sensor_data_sgmntd, sensation_data, IMU_map, M_sntn_Map,...
    ext_backward, ext_forward, FM_dilation_time, Fs_sensor, Fs_sensation)

% MATCH_WITH_M_SENSATION Summary of this function goes here
%   Input variables:  sensor_data_sgmntd- a cell variable with single cell/multiple cells
%                                         Each cell contains data from a sensor or a combination.
%                     sensation_data, IMU_map, M_sntn_Map- cell variables with single cell   
%                     ext_bakward, ext_forward, FM_dilation_time- scalar values 
%                     Fs_sensor, Fs_sensation- scalar values
%   Output variables: TPD, FPD, TND, FND- vectors with number of rows equal to the 
%                                         number of cells in the sensor_data_sgmntd


% Calculation of fixed values
n_sensors = length(sensor_data_sgmntd);
matching_window_size = ext_backward + ext_forward; % window size is equal to the window size used to create the maternal sensation map
minm_overlap_time = FM_dilation_time/2; % Minimum overlap in second
DLB = round(ext_backward*Fs_sensor); % backward extension length in sample number
DLF = round(ext_forward*Fs_sensor); % forward extension length in sample number

% Variable decleration
TPD = zeros(n_sensors,1); % True positive detection
FND = zeros(n_sensors,1); % False negative detection
TND = zeros(n_sensors,1); % True negative detection
FPD = zeros(n_sensors,1); % False positive detection

% Labeling sensation data and determining number of maternal sensation detection
sensation_data_labeled = bwlabel(sensation_data);
n_Maternal_detected_movement = length(unique(sensation_data_labeled)) - 1; % 1 is deducted to remove the initial value of 0

% Loop for matching individual sensors
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
            
            if (~X)  % true when there is no coincidence
                
                % TPD and FND calculation for indivudial sensors
                Y = sum(sensor_data_sgmntd{j}.*indv_sensation_map); % Non-zero value gives the matching
                if (Y) % true if there is a coincidence
                    TPD(j) = TPD(j) + 1; % TPD incremented
                else
                    FND(j) = FND(j) + 1; % FND incremented
                end
                               
            end
            
        end
        
    end
    %

    % ------------------- Determination of TND and FPD  ------------------%    
    % Removal of the TPD and FND parts from the individual sensor data
    sensor_data_sgmntd_labeled = bwlabel(sensor_data_sgmntd{j});
    curnt_matched_vector = sensor_data_sgmntd_labeled.*M_sntn_Map; % Non-zero elements gives the matching. In M-sntn_Map multiple windows can overlap, which was not the case in case of sensation_data 
    curnt_matched_label = unique(curnt_matched_vector); % Gives the label of the matched sensor data segments    
    Arb_val = 4; % An arbitrary value
    
    if (length(curnt_matched_label)>1)
        curnt_matched_label = curnt_matched_label(2:end); % Removes the first element, which is 0        
        for m = 1 : length(curnt_matched_label)
            sensor_data_sgmntd{j}(sensor_data_sgmntd_labeled == curnt_matched_label(m)) = Arb_val; 
            % Assigns an arbitrary value to the TPD segments of the segmented signal
        end        
    end
    
    sensor_data_sgmntd{j}(M_sntn_Map == 1) = Arb_val; % Assigns an arbitrary value to the area under the M_sntn_Map
    sensor_data_sgmntd_removed = sensor_data_sgmntd{j}(sensor_data_sgmntd{j} ~= Arb_val); % Removes all the elements with value = Arb_val from the segmented data
    %
    
    % Calculation of TND and FPD for individual sensors
    L_removed = length(sensor_data_sgmntd_removed);
    index_window_start = 1;
    index_window_end = min(index_window_start + Fs_sensor*matching_window_size, L_removed);
    %
    while (index_window_start < L_removed)
        indv_window = sensor_data_sgmntd_removed(index_window_start: index_window_end);
        index_non_zero = find(indv_window);
        
        if length(index_non_zero) >= (minm_overlap_time*Fs_sensor)
            FPD(j) = FPD(j) + 1;
        else
            TND(j) = TND(j) + 1;
        end
        
        index_window_start = index_window_end + 1;
        index_window_end = min(index_window_start + Fs_sensor*matching_window_size, L_removed);
    end
 
end
%

end

