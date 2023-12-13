clc;
clear;
close all;

avg_N = 1000;
sigma_arr = [0.05 0.02 0.01];
tx_power_arr = [round(power(10,0),5) round(power(10,0.5),5) round(power(10,1),5) round(power(10,1.5),5) round(power(10,2),5)];
avg_mse_robust = zeros(3,5);
iteration = 1000;

% Then construct a ParforProgressbar object:
% ppm = ParforProgressbar(iteration);
% ppm = ParforProgressbar(iteration, 'showWorkerProgress', true);


for i=1:avg_N
    temp = zeros(3,5);
    for sig = 1:length(sigma_arr)
        for power_index = 1:length(tx_power_arr)
            [mse,~]=mmse(40,sigma_arr(sig),tx_power_arr(power_index),1);
            temp(sig,power_index) = mse;
        end
    end
    avg_mse_robust = avg_mse_robust+temp;

%     pause(100/iteration);
%     % increment counter to track progress
%     ppm.increment();
end

% Delete the progress handle when the parfor loop is done (otherwise the timer that keeps updating the progress might not stop).
% delete(ppm);



avg_mse_robust = avg_mse_robust./avg_N;

