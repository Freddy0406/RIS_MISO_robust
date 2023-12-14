clc;
clear;
close all;

avg_N = 1000;                         % Average times
variance_arr = [0.05 0.02 0.01];    % Figure 1 with different variance
tx_power_arr = [round(power(10,0),5) round(power(10,0.5),5) round(power(10,1),5) round(power(10,1.5),5) round(power(10,2),5)];
% Figure 1 with different transmit power

avg_mse_robust = zeros(3,5);
avg_mse_nonrobust = zeros(3,5);
iteration = avg_N;                  %parfor iteration


%% To make figure 1

ppm = ParforProgressbar(iteration);
ppm = ParforProgressbar(iteration, 'showWorkerProgress', true);
parfor i=1:avg_N
    temp_robust = zeros(3,5);
    temp_nonrobust = zeros(3,5);
    for sig = 1:length(variance_arr)
        for power_index = 1:length(tx_power_arr)
            [mse,~]=mmse(40,variance_arr(sig),tx_power_arr(power_index),1,0);
            temp_robust(sig,power_index) = mse;
            [mse,~]=mmse(40,variance_arr(sig),tx_power_arr(power_index),2,0);
            temp_nonrobust(sig,power_index) = mse;
        end
    end
    avg_mse_robust = avg_mse_robust+temp_robust;
    avg_mse_nonrobust = avg_mse_nonrobust+temp_nonrobust;
    pause(100/iteration);
    ppm.increment();
end
delete(ppm);

avg_mse_robust = avg_mse_robust./avg_N;
avg_mse_nonrobust = avg_mse_nonrobust./avg_N;

