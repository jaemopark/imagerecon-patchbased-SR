load sample_data.mat
%% use number of iteration = 50 as threshold to end the algorithm
% target_th_type = 0,target_th = 50
[c_super_res,e_final,iteration_final] = super_res_2d(brain_13C_lac,brain_1H,possibility,5,64,...
                                   0,5);
%% ues reconstruction error < 1 as threshold to end the algorithm
% target_th_type = 1,target_th = 50
[c_super_res,e_final,iteration_final] = super_res_2d(brain_13C_lac,brain_1H,possibility,5,64,...
                                   1,1);