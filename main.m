clc;
clear;
close all;

avg_N = 100;
sigma_arr = [0.05 0.02 0.01];
avg_mse = zeros(2,length(sigma_arr));

for a = 1:length(sigma_arr)
    for b = 1:avg_N
        initial(40,sigma_arr(a),power(10,0.2));
        [mse,~]=mmse(1);
        avg_mse(1,a) = avg_mse(1,a)+mse;
        [mse,~]=mmse(2);
        avg_mse(2,a) = avg_mse(2,a)+mse;
    end
end

avg_mse = avg_mse./avg_N;

