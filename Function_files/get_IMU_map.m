function [IMU_map] = get_IMU_map(IMU_data, data_file_name, Fs_sensor)

% Summary of this function goes here
%   Input variables:  IMU_data- A column vector of IMU data 
%                     data_file_names- a char variable with data file name
%   Output variables: IMU_map- A column vector with segmented IMU data 


% Counting the number of data files to be loaded
%   data_file_names is a string in case of a single data file selected.
%   Otherwise, it is a cell variable containing 1 namein each cell    

IMU_threshold = [0.003 0.002]; % fixed threshold value obtained through seperate testing
IMU_dilation_time = 4.0; % dilation length in seconds
IMU_dilation_size = round(IMU_dilation_time*Fs_sensor); % dilation lenght in sample number
SE_IMU = strel('line', IMU_dilation_size, 90); % Creates a linear structuring element (vertical, as deg=90) that will have values 1 with a length of dilation_size;

% --------- Segmentaiton of IMU data and creation of IMU_map ---------%
if data_file_name(2) == '3'
    IMU_map = abs(IMU_data) >= IMU_threshold(2); % Threshold for Subject 3
else
    IMU_map = abs(IMU_data) >= IMU_threshold(1); % Threshold for Subject 1 or 2
end

% ----------------------- Dilation of IMU data ------------------------
IMU_map = imdilate(IMU_map, SE_IMU); % Dilate or expand the ROI's(points with value = 1) by dilation_size (half above and half below),as defined by SE
%

end

