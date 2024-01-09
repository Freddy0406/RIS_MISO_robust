clc;
clear;
close all;

t=1;                    %Set times
avg_N = 1000;                         % Average times
variance_arr = [0.05 0.01 0.002];    % Figure 1 with different variance
tx_power_arr = [round(power(10,0),5) round(power(10,0.5),5) round(power(10,1),5) round(power(10,1.5),5) round(power(10,2),5)];
% Figure 1 with different transmit power
%% Initialization parameter
N = 40;
P0 = tx_power_arr(1);
variance = variance_arr(1);
M = 4;                          %BS antenna amount
L0 = power(10,(-30/10));        %Path loss -30dB => linear 
AP_loc = [0 0];                 %AP location
IRS_loc = [100,0];              %IRS location
UE_loc = [100 20];              %User location

%% Initialize s
s = randn(1,1,"like",1i);
c = zeros(1,100);

AP_UE_dis = sqrt(sum((UE_loc-AP_loc).^2));              %AP-UE distance  
AP_IRS_dis = sqrt(sum((AP_loc-IRS_loc).^2));            %AP-IRS distance
UE_IRS_dis = sqrt(sum((UE_loc-IRS_loc).^2));            %UE-IRS distance

alpha_LoS = -2;                                         %Path Loss Exponent of LoS(AP->IRS->UE)
alpha_nLoS = -3;                                        %Path Loss Exponent of non-LoS(AP->UE)
LoS_PL =  L0*power((AP_IRS_dis+UE_IRS_dis),alpha_LoS);  %Path Loss of LoS(AP->IRS->UE)
nLoS_PL = L0*power((AP_UE_dis),alpha_nLoS);             %Path Loss of non-LoS(AP->UE)
epsilon = power(10,-4);                                 %Limit of optimization iteration



%Fading Channel
Rayleigh_fading = sqrt(0.5) * (randn(size(s)) + 1i * randn(size(s)));       %Random Rayleigh fading noise
                                                                            %Use in non-LoS(AP->UE)

K = 10;                                                                     %Rician fading factor
LoS = sqrt(K / (K + 1));                                                    %Rician fading main lobe
nLoS = sqrt(1 / (K + 1)) * (randn(size(s)) + 1i * randn(size(s)));          %Rician fading side lobe

%Noise
noise_var = power(10,-11);                      %var_n = 110dBm = 10^-11 mW
noise = sqrt(noise_var)*randn(1,1,"like",1i);   %Generate random 1X1 noise

% Initialize w 
w = zeros(M,1,100);
for i = 1:M
    w(i,1,t) = sqrt(P0)/sqrt(M);
end

% Initialize theta (Start with random phases)
btheta = zeros(N,N,100);
min = 0;
max = 2*pi;
v = exp(1i*(min+rand(N,1)*(max-min)));
btheta(:,:,t) = diag(v);

% Generate delta_G, delta_hd, delta_hr (AP-IRS, IRS-User, AP-User link)
delta_G = sqrt(variance)*randn(N,M,"like",1i);
delta_hd = sqrt(variance)*randn(M,1,"like",1i);
delta_hr = sqrt(variance)*randn(N,1,"like",1i);

% Emulate estimation channel

G = vertcat(eye(M),zeros((N-M),M));     %Ideal channel => Identity matrix
hd = randn(M,1)+1i*randn(M,1);          %Ideal channel => Rayleigh fading                         
hr = ones(N,1);                                   %Ideal channel => All ones

G_hat = G-delta_G;          %G = G_hat+delta_G ; G_hat=G-delta_G        (AP-IRS)
hd_hat = hd-delta_hd;       %hd = hd_hat+delta_hd ; hd_hat=hd-delta_hd  (IRS-User)
hr_hat = hr-delta_hr;       %hr = hr_hat+delta_hr ; hr_hat=hr-delta_hr  (AP-User)



%% Start optimization   
mse_before =0.0;


while true
    fprintf('The %d generation...\n',t);
    %% Optimization of c
    fprintf('\tOptimization of c...\n');        
    A = (G_hat'*btheta(:,:,t)'*hr_hat)*hr_hat'*btheta(:,:,t)*G_hat+...          %According to the formula              
        hd_hat*hr_hat'*btheta(:,:,t)*G_hat+G_hat'*btheta(:,:,t)'*hr_hat*hd_hat'+...
        hd_hat*hd_hat'+variance*power(norm(hr_hat),2)*eye(M)+(variance*G_hat')*G_hat+...
        (N*variance*variance+variance)*eye(M);

    alpha = (G_hat'*btheta(:,:,t)'*hr_hat+hd_hat);

    c(t) = (w(:,:,t)'*alpha)/(w(:,:,t)'*A*w(:,:,t)+noise_var);          %Wiener Filter

    %% Optimization of w
    fprintf('\tOptimization of w...\n');
    lambda=0;
    w(:,:,t+1) = (power(abs(c(t)),2)*A+lambda*eye(M))^-1*(alpha*conj(c(t)));  %Assume lambda(Lagrange mult.) = 0
    if(power(norm(w(:,:,t+1)),2)>P0)
        [lambda] = search(c(t),alpha,A,P0,M);
        w(:,:,t+1) = (power(abs(c(t)),2)*A+lambda*eye(M))^-1*(alpha*conj(c(t)));
    else
        lambda=0;
        w(:,:,t+1) = (power(abs(c(t)),2)*A+lambda*eye(M))^-1*(alpha*conj(c(t)));
    end
    
    %% Optimization of theta
    fprintf('\tOptimization of btheta...\n');
    Phi = diag(hr_hat')*G_hat*w(:,:,t+1)*c(t);
    d = hd_hat'*w(:,:,t+1)*c(t);
    Q = Phi*Phi';
    q = Phi*(1-conj(d)');
    v = diag(btheta(:,:,t));
    [new_v] = opt_v(Q,q,v);
    v = new_v;
    btheta(:,:,t+1) = diag(v);

    %% Calculate mse
    fprintf('\tCalculate mmse...\n\n');
    btheta_2D = btheta(:,:,t+1);
    mse = power(abs(c(t)),2)*(w(:,:,t+1)'*A*w(:,:,t+1)+noise_var)-w(:,:,t+1)'*alpha*conj(c(t))-c(t).*alpha'*w(:,:,t+1)+1;

    if((mse_before-mse)>epsilon)
        mse_before = mse;
    elseif((mse_before-mse)<=epsilon)
        break;
    end              
    t = t+1;
end

