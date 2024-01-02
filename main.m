clc;
clear;
close all;


avg_N = 1000;                           % Average times
variance_arr = [0.05 0.01 0.002];       % Figure 1 with different variance
tx_power_arr = [round(power(10,0),5) round(power(10,0.5),5) round(power(10,1),5) round(power(10,1.5),5) round(power(10,2),5)];
% Figure 2 with different transmit power

avg_mse_robust = zeros(3,5);
avg_mse_nonrobust = zeros(3,5);
iteration = 1000;                       %parfor iteration

% Store data for figure 3
variance_arr_2 = [0.08 0.04];
N_array = [10 20 30 40 50 60];
avg_mse_robust_2 = zeros(2,6);
avg_mse_nonrobust_2 = zeros(2,6);
avg_mse_robust_discrete = zeros(2,6);
avg_mse_without_IRS = zeros(1,6);



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
            [mse,~]=mmse(40,variance_arr(sig),tx_power_arr(power_index),2,0);
            temp_nonrobust(sig,power_index) = mse;
        end
    end
    avg_mse_robust = avg_mse_robust+temp_robust;
    avg_mse_nonrobust = avg_mse_nonrobust+temp_nonrobust;
    

%     temp_mse_robust_2 = zeros(2,6);
%     temp_mse_nonrobust_2 = zeros(2,6);
%     temp_mse_robust_discrete = zeros(2,6);
%     temp_mse_without_IRS = zeros(1,6);
% 
%     
%     for sig = 1:length(variance_arr_2)
%         for N = 1:length(N_array)
%             [mse,~]=mmse(N_array(N),variance_arr_2(sig),round(power(10,1),5),1,0);
%             temp_mse_robust_2(sig,N) = mse;
%             [mse,~]=mmse(N_array(N),variance_arr_2(sig),round(power(10,1),5),2,0);
%             temp_mse_nonrobust_2(sig,N) = mse;
%         end
%     end
%     avg_mse_robust_2 = avg_mse_robust_2+temp_mse_robust_2;
%     avg_mse_nonrobust_2 = avg_mse_nonrobust_2+temp_mse_nonrobust_2;
% 
% 
%     for N = 1:length(N_array)
%         [mse,~]=mmse(N_array(N),variance_arr_2(1),round(power(10,1),5),3,2);
%         temp_mse_robust_discrete(1,N) = mse;
%         [mse,~]=mmse(N_array(N),variance_arr_2(1),round(power(10,1),5),3,3);
%         temp_mse_robust_discrete(2,N) = mse;
%         [mse,~]=mmse(N_array(N),variance_arr_2(1),round(power(10,1),5),4,0);
%         temp_mse_without_IRS(1,N) = mse;
%     end
% 
%     avg_mse_robust_discrete = avg_mse_robust_discrete+temp_mse_robust_discrete;
%     avg_mse_without_IRS = avg_mse_without_IRS+temp_mse_without_IRS;

    %Progressbar
    pause(100/iteration);
    ppm.increment();
end

%Delete Progressbar
delete(ppm);

%% Plot for figure2

%Average result
avg_mse_robust = avg_mse_robust./avg_N;
avg_mse_nonrobust = avg_mse_nonrobust./avg_N;

%Change into dBm
tx_power_arr = 10*log10(tx_power_arr);

figure(1)
semilogy(tx_power_arr,avg_mse_robust(1,:),'r-o',tx_power_arr,avg_mse_robust(2,:),'b-square',tx_power_arr,avg_mse_robust(3,:),'k-diamond',tx_power_arr,avg_mse_nonrobust(1,:),'r--x',tx_power_arr,avg_mse_nonrobust(2,:),'g--^',tx_power_arr,avg_mse_nonrobust(3,:),'m--+');
xlabel("Transmit Power (dBm)")
ylabel("MSE")
legend('Robust, {\sigma^2}=0.05','Robust, {\sigma^2}=0.01','Robust, {\sigma^2}=0.002','Non-Robust, {\sigma^2}=0.05','Non-Robust, {\sigma^2}=0.01','Non-Robust, {\sigma^2}=0.002')

%% Plot for figure3

%Average result
avg_mse_robust_2 = avg_mse_robust_2./avg_N;
avg_mse_nonrobust_2 = avg_mse_nonrobust_2./avg_N;
avg_mse_robust_discrete = avg_mse_robust_discrete./avg_N;
avg_mse_without_IRS = avg_mse_without_IRS./avg_N;

figure(2);
semilogy(N_array,avg_mse_robust_2(1,:),'r-o',N_array,avg_mse_robust_2(2,:),'m-square',N_array,avg_mse_nonrobust_2(1,:),'c-*',N_array,avg_mse_nonrobust_2(2,:),'g-|',N_array,avg_mse_robust_discrete(1,:),'k--diamond',N_array,avg_mse_robust_discrete(2,:),'b--^',N_array,avg_mse_without_IRS,'k--x');
xlabel("Reflecting Elements")
ylabel("MSE")
legend('Robust, {\sigma^2}=0.08','Robust, {\sigma^2}=0.04','Non-Robust, {\sigma^2}=0.08','Non-Robust, {\sigma^2}=0.04','2-bit IRS, {\sigma^2}=0.08','3-bit IRS, {\sigma^2}=0.08','without IRS, {\sigma^2}=0.08')


