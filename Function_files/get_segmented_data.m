function [sensor_data_sgmntd, h] = get_segmented_data(sensor_data, min_SN, IMU_map, dilation_time, Fs_sensor)

% GET_SEGMENTED_DATA Summary of this function goes here
%   Input variables:  sensor_data- a cell variable;
%                     min_SN- a vector/ a scalar;
%                     IMU_map- a column vector
%                     Fs_sensor, FM_dilation_time- a scalar
%   Output variables: sensor_data_sgmntd- a cell variable of same size as the input variable sensor_data_fltd;
%                     h- a vector
%

low_signal_quantile = 0.25;
sgmntd_signal_cutoff = [0.0001 0.0001 0.0001 0.0001 0.0001 0.0001];
dilation_size = round (dilation_time * Fs_sensor); 
SE = strel('line', dilation_size, 90); % linear element necessary for dilation operation

n_sensors = length(sensor_data);

if(isscalar(min_SN))
    min_SN_new (1 : n_sensors) = min_SN; % In the case where only a single value is given for FM_min_SN, it is expanded for all the sensors
else
    min_SN_new = min_SN;
end

h = zeros(1,n_sensors); % Variable for thresold
sensor_data_sgmntd = cell(1,n_sensors);

% Thresholding
for j = 1 : n_sensors

    % Determining the threshold
    s = abs(sensor_data{j});
    LQ = quantile(s,low_signal_quantile); % Rerurns the quantile value for low 25% (= low_signal_quantile) of the signal
    e = s(s <= LQ); % Signal noise
    h(j) = min_SN_new(j)*median(e); % Threshold value. Each row will contain threshold value for each data file
    
    if isnan(h(j)) % Check if h = Nan. This happens when e=[], as median(e)= NAN for that case!!!
        h(j) = Inf;
    end    
    if h(j) < sgmntd_signal_cutoff(j)
        h(j) = Inf; % Precaution against too noisy signal
    end

    % Thresholding    
    sensor_data_sgmntd{j} = (s >= h(j)); % Considering the signals that are above the threshhold value; data point for noise signal becomes zero and singal become 1

    % Exclusion of body movement data
    sensor_data_sgmntd{j} = sensor_data_sgmntd{j}.*(1-IMU_map); % Exclusion of body movement data

    % Dilation of the thresholded data
    sensor_data_sgmntd{j} = imdilate(sensor_data_sgmntd{j}, SE); % For individual sensor performance

end

end

