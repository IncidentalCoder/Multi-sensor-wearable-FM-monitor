function [SEN_all,PPV_all,SPE_all,ACC_all,FS_all,FPR_all] = get_performance_params(TPD_all,FPD_all,TND_all,FND_all)

% GET_PERFORMANCE_PARAMS Summary of this function goes here
%   Input variables:  TPD_all, FPD_all, TND_all, FND_all- single cell/multi-cell variable.
%                     Number of cell indicates number of sensor data or
%                     combination data provided together.
%                     Each cell containes a vector with row size = n_data_files 
%
%   Output variables: SEN_all,PPV_all,SPE_all,ACC_all,FS_all,FPR_all- cell variable with
%                     size same as the input variables.        

n_sensors = length(TPD_all); % TPD_all, FPD_all, TND_all, and FND_all all are will be of same length

SEN_all = cell(1,n_sensors);
PPV_all = cell(1,n_sensors);
SPE_all = cell(1,n_sensors);
ACC_all = cell(1,n_sensors);
FS_all = cell(1,n_sensors);
FPR_all = cell(1,n_sensors);

B = 1; % Beta value for F_B score calculation.

for i = 1 : n_sensors
    
    SEN_all{i} = TPD_all{i}./(TPD_all{i} + FND_all{i}); % Operates on individual data sets
    PPV_all{i} = TPD_all{i}./(TPD_all{i} + FPD_all{i});
    SPE_all{i} = TND_all{i}./(TND_all{i} + FPD_all{i});
    ACC_all{i} = (TPD_all{i} + TND_all{i})./(TPD_all{i} + FND_all{i} + FPD_all{i} + TND_all{i});
    FS_all{i} = (1+B^2) * (PPV_all{i} .* SEN_all{i})./(B^2*PPV_all{i} + SEN_all{i});
    
    FPR_all{i} = FPD_all{i}./(FPD_all{i} + TND_all{i});
    
end


end

