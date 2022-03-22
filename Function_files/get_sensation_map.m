function [M_sntn_map] = get_sensation_map(sensation_data, IMU_map, ext_backward, ext_forward, Fs_sensor, Fs_sensation)

% Summary of this function goes here
%   Input variables-  sensation_data- a column vector with the sensation data
%                     IMU_map- a column vector
%                     ext_backward, ext_forward- scalar values
%                     Fs_sensor,Fs_sensation- scalar values    %                     
%   Output variables: M_sntn_map- a column vector with all the sensation windows joined together  

M_event_index = find(sensation_data); % Sample numbers for maternal sensation detection
M_sntn_map = zeros(length(IMU_map),1); % Initializing the map with zeros everywhere

% Parameters for creating M_sensation_map
DLB = round(ext_backward*Fs_sensor); % backward extension length in number of sample
DLF = round(ext_forward*Fs_sensor); % forward extension length in number of sample

for j = 1:length(M_event_index) % M_event_index contains index of all the maternal sensation detections

    % Getting the index values for the map
    L = M_event_index(j); % Sample no. corresponding to a maternal sensation
    L1 = L*round(Fs_sensor/Fs_sensation) - DLB; % sample no. for the starting point of this sensation in the map
    L2 = L*round(Fs_sensor/Fs_sensation) + DLF; % sample no. for the ending point of this sensation in the map
    L1 = max(L1, 1); % Just a check so that L1 remains higher than 1st data sample
    L2 = min(L2, length(M_sntn_map)); % Just a check so that L2 remains lower than the last data sample
    %
    % Generating the map- a single vector with all the sensation data mapping
    M_sntn_map(L1:L2) = 1; % Assigns 1 to the location defined by L1:L2

    % Removal of the maternal sensation that has coincided with the body movement
    X = sum(M_sntn_map(L1:L2).*IMU_map(L1:L2)); % this is non-zeor if there is a coincidence
    if (X)
        M_sntn_map(L1:L2) = 0; % Removes the sensation data from the map
    end

end

end

