clc;
clear;
close all;


avg_N = 1000;                           % Average times
variance_arr = [0.05 0.01 0.002];       % Figure 1 with different variance
tx_power_arr = [round(power(10,0),5) round(power(10,0.5),5) round(power(10,1),5) round(power(10,1.5),5) round(power(10,2),5)];
% Figure 2 with different transmit power

avg_mse_robust = zeros(3,5);
iteration = 1000;                       %parfor iteration




%% Simulation start

%Progressbar
ppm = ParforProgressbar(iteration);
ppm = ParforProgressbar(iteration, 'showWorkerProgress', true);


parfor i=1:1000

    % Data for figure2
    temp_robust = zeros(3,5);
    temp_nonrobust = zeros(3,5);
    for sig = 1:length(variance_arr)
        for power_index = 1:length(tx_power_arr)
            [mse,~]=mmse(40,variance_arr(sig),tx_power_arr(power_index),1,0);
            temp_robust(sig,power_index) = mse;
        end
    end

    %Progressbar
    pause(100/iteration);
    ppm.increment();
end

%Delete Progressbar
delete(ppm);

%% Plot for figure2

%Average result
avg_mse_robust = avg_mse_robust./avg_N;

%Change into dBm
tx_power_arr = 10*log10(tx_power_arr);

figure(2)
semilogy(tx_power_arr,avg_mse_robust(1,:),'r-o',tx_power_arr,avg_mse_robust(2,:),'b-square',tx_power_arr,avg_mse_robust(3,:));
xlabel("Transmit Power (dBm)")
ylabel("MSE")
legend('Robust, {\sigma^2}=0.05','Robust, {\sigma^2}=0.01','Robust, {\sigma^2}=0.002')

