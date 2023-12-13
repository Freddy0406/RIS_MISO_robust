function initial(N,variance,P0)
    %% Initialization parameter
    %Lr         %Reflecting element amount per RIS
    %PT         %Maximum transmit power of BS
    %var        %var_g=var_r=var_d=var

    M = 4;                          %BS antenna amount
    s = randn(1,1,"like",1i);       %Symbol for user(zero-meam, unit power=>variance=1)
    L0 = -30;                           
    AP_loc = [0 0];                 %AP location
    IRS_loc = [100,0];              %IRS location
    UE_loc = [100 20];              %User location
    AP_UE_dis = sqrt(sum((UE_loc-AP_loc).^2));
    AP_IRS_dis = sqrt(sum((AP_loc-IRS_loc).^2));
    UE_IRS_dis = sqrt(sum((UE_loc-IRS_loc).^2));
    alpha_LoS = -2;
    alpha_nLoS = -3;
    LoS_PLE =  L0*power((AP_IRS_dis+UE_IRS_dis),alpha_LoS);
    nLoS_PLE = L0*power((AP_UE_dis),alpha_nLoS);
    epsilon = power(10,-4);

    noise_var = power(10,-11);
    noise = sqrt(noise_var)*randn(1,1,"like",1i);


    %% Initialize w
    w = zeros(M,1);

    for i = 1:M
        w(i) = sqrt(P0)/sqrt(M);
    end

    %% Initialize theta
    min = 0;
    max = 2*pi;
    v = exp(1i*(min+rand(N,1)*(max-min)));
    btheta = diag(v);

    %% delta_G, delta_hd, delta_hr (AP-IRS, IRS-User, AP-User link)
    
    delta_G = sqrt(variance)*randn(N,M,"like",1i);
    delta_hd = sqrt(variance)*randn(M,1,"like",1i);
    delta_hr = sqrt(variance)*randn(N,1,"like",1i);
    
    %% Emulate estimation channel

    G = vertcat(eye(M),zeros((N-M),M));      %Ideal channel
    hd = ones(M,1);                         
    hr = ones(N,1);
    
    G_hat = G-delta_G;
    hd_hat = hd-delta_hd;
    hr_hat = hr-delta_hr;

    save('workspace');
end