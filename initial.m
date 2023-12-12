function initialization(N,variance,P0)
    %% Initialization parameter
    %Lr         %Reflecting element amount per RIS
    %PT         %Maximum transmit power of BS
    %var        %var_g=var_r=var_d=var

    M = 4;                          %BS antenna amount
    s = randn(1,1,"like",1i);       %Symbol for user(zero-meam, unit power=>variance=1)
    AP_loc = [0 0];                 %AP location
    IRS_loc = [100,0];              %IRS location
    UE_loc = [100 20];              %User location
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

    G = sqrt(variance)*randn(N,M,"like",1i);
    hd = sqrt(variance)*randn(M,1,"like",1i);
    hr = sqrt(variance)*randn(N,1,"like",1i);
    
    G_hat = G-delta_G;
    hd_hat = hd-delta_hd;
    hr_hat = hr-delta_hr;

    save('workspace');
end