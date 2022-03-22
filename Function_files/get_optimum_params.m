function [r_best, FM_min_SN_best, Threshold_best, TPD_best, FPD_best, TND_best, FND_best, SEN_best, PPV_best, SPE_best, ACC_best, FS_best, FPR_best] ...
    = get_optimum_params(n_data_files, n_sensors, FM_min_SN_initial, SN_increment, TPD_all, FPD_all, TND_all, FND_all, SEN_all, PPV_all, SPE_all, ACC_all, ...
    FS_all, FPR_all, h)

%GET_OPTIMUM_PARAMS Summary of this function goes here
%   Detailed explanation goes here

% Variable decleration
r_best = zeros(n_data_files, n_sensors); % stores the iteration number that has the maximum value of F-score
FM_min_SN_best = zeros(n_data_files, n_sensors); % Stores the best SN ratio for each sensor in each data set
Threshold_best = zeros(n_data_files, n_sensors); % Stores the best threshold value for each sensor in each data file

TPD_best = zeros(n_data_files, n_sensors);
FPD_best = zeros(n_data_files, n_sensors);
TND_best = zeros(n_data_files, n_sensors);
FND_best = zeros(n_data_files, n_sensors);
FS_best = zeros(n_data_files, n_sensors);
SEN_best = zeros(n_data_files, n_sensors);
PPV_best = zeros(n_data_files, n_sensors);
SPE_best = zeros(n_data_files, n_sensors);
ACC_best = zeros(n_data_files, n_sensors);
FPR_best = zeros(n_data_files, n_sensors);

% for i = 1 : n_sensors
%     
%     [~, r_best(:,i)] = max(FS_all{i}, [], 2); % iteration number with best FS- each data set along each row and each sensor along each column
%     FM_min_SN_best(:,i) = FM_min_SN_initial + (r_best(:, i) - 1) * SN_increment;    
%     
% end

for i = 1 : n_data_files

    for j = 1 : n_sensors

        [~, r_best(i,j)] = max(FS_all{j}(i,:), [], 2); % Fs_all{j} gives FS for j-th sensors. Then the max gives the max of i-th row
        %         [~, r_best(i,j)] = max(ACC_all{j}(i,:), [], 2); % To get the optimization for Accuracy or PABAK

        FM_min_SN_best(i,j) = FM_min_SN_initial + (r_best(i,j) - 1) * SN_increment;

        TPD_best(i,j) = TPD_all{j}(i,r_best(i,j)); % Each data set along each row and each sensor along each column
        FPD_best(i,j) = FPD_all{j}(i,r_best(i,j));
        TND_best(i,j) = TND_all{j}(i,r_best(i,j));
        FND_best(i,j) = FND_all{j}(i,r_best(i,j));

        SEN_best(i,j) = SEN_all{j}(i,r_best(i,j));
        PPV_best(i,j) = PPV_all{j}(i,r_best(i,j));
        SPE_best(i,j) = SPE_all{j}(i,r_best(i,j));
        ACC_best(i,j) = ACC_all{j}(i,r_best(i,j));
        FS_best (i,j) = FS_all{j}(i,r_best(i,j));
        FPR_best (i,j) = FPR_all{j}(i,r_best(i,j));
        
        Threshold_best(i,j) = h{r_best(i,j)}(i,j); % in h, each iteration is saved as a cell file. r_best(i,j) gives the best iteration number 
        
    end
    
end
%

% for i = 1:n_FM_sensors
%     
%     % For individual sensors
%     [FS_best_overall(i), r_best_overall(i)] = max(FS_all_overall{i});   
%     FM_min_SN_best_overall(i) = FM_min_SN_initial + (r_best_overall(i)-1) * SN_increment;
%     
%     TPD_best_overall (i) = TPD_all_overall{i}(r_best_overall(i));
%     FPD_best_overall (i) = FPD_all_overall{i}(r_best_overall(i));
%     TND_best_overall (i) = TND_all_overall{i}(r_best_overall(i));
%     FND_best_overall (i) = FND_all_overall{i}(r_best_overall(i));
%     SEN_best_overall (i) = SEN_all_overall{i}(r_best_overall(i));
%     PPV_best_overall (i) = PPV_all_overall{i}(r_best_overall(i));
%     SPE_best_overall (i) = SPE_all_overall{i}(r_best_overall(i));
%     ACC_best_overall (i) = ACC_all_overall{i}(r_best_overall(i));
%     
%      % For combined sensors
%     if i <= n_FM_sensors/2       
%         [FS_best_cmbd_overall(i), r_best_cmbd_overall(i)] = max(FS_all_cmbd_overall{i});
%         FM_min_SN_best_cmbd_overall(i) = FM_min_SN_initial + (r_best_cmbd_overall(i)-1) * SN_increment;
%         
%         TPD_best_cmbd_overall (i) = TPD_all_cmbd_overall{i}(r_best_cmbd_overall(i));
%         FPD_best_cmbd_overall (i) = FPD_all_cmbd_overall{i}(r_best_cmbd_overall(i));
%         TND_best_cmbd_overall (i) = TND_all_cmbd_overall{i}(r_best_cmbd_overall(i));
%         FND_best_cmbd_overall (i) = FND_all_cmbd_overall{i}(r_best_cmbd_overall(i));  
%         SEN_best_cmbd_overall (i) = SEN_all_cmbd_overall{i}(r_best_cmbd_overall(i));
%         PPV_best_cmbd_overall (i) = PPV_all_cmbd_overall{i}(r_best_cmbd_overall(i));
%         SPE_best_cmbd_overall (i) = SPE_all_cmbd_overall{i}(r_best_cmbd_overall(i));
%         ACC_best_cmbd_overall (i) = ACC_all_cmbd_overall{i}(r_best_cmbd_overall(i));      
%     end
    
% end


end

