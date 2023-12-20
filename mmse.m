function [mse,t]=mmse(N,variance,P0,mode,bit_of_phase)
%% Initialization parameter
    %N  Reflecting element amount per RIS
    %P0 Maximum transmit power of BS
    %variance var_g=var_r=var_d=var
    M = 4;                          %BS antenna amount
    L0 = -30;                       %Path loss dB -30dB         
    AP_loc = [0 0];                 %AP location
    IRS_loc = [100,0];              %IRS location
    UE_loc = [100 20];              %User location

    s = rand;
    if(s>0.5)
        s = 1;
    else 
        s = -1;
    end

    AP_UE_dis = sqrt(sum((UE_loc-AP_loc).^2));              %AP-UE distance  
    AP_IRS_dis = sqrt(sum((AP_loc-IRS_loc).^2));            %AP-IRS distance
    UE_IRS_dis = sqrt(sum((UE_loc-IRS_loc).^2));            %UE-IRS distance
    alpha_LoS = -2;                                         %Path Loss Exponent of LoS(AP->IRS->UE)
    alpha_nLoS = -3;                                        %Path Loss Exponent of non-LoS(AP->UE)
    LoS_PL =  L0*power((AP_IRS_dis+UE_IRS_dis),alpha_LoS);  %Path Loss of LoS(AP->IRS->UE)
    LoS_PL = 10^(LoS_PL/10); 
    nLoS_PL = L0*power((AP_UE_dis),alpha_nLoS);             %Path Loss of non-LoS(AP->UE)
    nLoS_PL = 10^(nLoS_PL/10);
    epsilon = power(10,-4);                                 %Limit of optimization iteration


    %Fading Channel
    Rayleigh_fading = sqrt(0.5) * (randn(size(s)) + 1i * randn(size(s)));       %Random Rayleigh fading noise
                                                                                %Use in non-LoS(AP->UE)

    K = 10;                                                                     %Rician fading factor
    K_linear = 10^(K/10);                                                       %Use in LoS(AP->IRS->UE)
    LoS = sqrt(K_linear / (K_linear + 1));                                      %Rician fading main lobe
    nLoS = sqrt(1 / (K_linear + 1)) * (randn(size(s)) + 1i * randn(size(s)));   %Rician fading side lobe

    %Noise
    noise_var = power(10,-11);                      %var_n = 110dBm = 10^-11 mW
    noise = sqrt(noise_var)*randn(1,1,"like",1i);   %Generate random 1X1 noise

    % Initialize w 
    w = zeros(M,1);
    for i = 1:M
        w(i) = sqrt(P0)/sqrt(M);
    end

    % Initialize theta (Start with random phases)
    min = 0;
    max = 2*pi;
    v = exp(1i*(min+rand(N,1)*(max-min)));
    btheta = diag(v);

    % Generate delta_G, delta_hd, delta_hr (AP-IRS, IRS-User, AP-User link)
    delta_G = sqrt(variance)*randn(N,M,"like",1i);
    delta_hd = sqrt(variance)*randn(M,1,"like",1i);
    delta_hr = sqrt(variance)*randn(N,1,"like",1i);
    
    % Emulate estimation channel

    G = vertcat(eye(M),zeros((N-M),M));     %Ideal channel => Identity matrix
    hd = ones(M,1);                         %Ideal channel => All ones                         
    hr = ones(N,1);                         %Ideal channel => All ones
    
    G_hat = G-delta_G;          %G = G_hat+delta_G ; G_hat=G-delta_G        (AP-IRS)
    hd_hat = hd-delta_hd;       %hd = hd_hat+delta_hd ; hd_hat=hd-delta_hd  (IRS-User)
    hr_hat = hr-delta_hr;       %hr = hr_hat+delta_hr ; hr_hat=hr-delta_hr  (AP-User)


    
    %% Start optimization   
    t=1;                    %Set times
%% mode 1:The proposed robust design
    if (mode==1)           
        while true
            fprintf('The %d generation...\n',t);
            %% Optimization of c
            fprintf('\tOptimization of c...\n');        
            A = (G_hat'*btheta'*hr_hat)*hr_hat'*btheta*G_hat+...          %According to the formula              
                hd_hat*hr_hat'*btheta*G_hat+G_hat'*btheta'*hr_hat*hd_hat'+...
                hd_hat*hd_hat'+variance*power(norm(hr_hat),2)*eye(M)+(variance*G_hat')*G_hat+...
                (N*variance*variance+variance)*eye(M);

            alpha = (G_hat'*btheta'*hr_hat+hd_hat);
        
            c = (w'*alpha)/(w'*A*w+noise_var);          %Wiener Filter
        
            %% Optimization of w
            fprintf('\tOptimization of w...\n');
            w =  (power(abs(c),2)*A)\(alpha*conj(c));  %Assume lambda(Lagrange mult.) = 0
            power(norm(w),2);
            
            if(power(norm(w),2)>P0)
                [lambda] = search(c,alpha,A,P0,M);
                w = (power(abs(c),2)*A+lambda*eye(M))\(alpha*conj(c));
            else
                lambda = 0;
                w = (power(abs(c),2)*A+lambda*eye(M))\(alpha*conj(c));
            end
            
            %% Optimization of theta
            fprintf('\tOptimization of btheta...\n');
            Phi = diag(hr_hat')*G_hat*w*c;
            d = hd_hat'*w*c;
            Q = Phi*Phi';
            q = Phi*(1-conj(d)');
            [new_v] = opt_v(Q,q,v);
            v = new_v;
            btheta = diag(v);
        
            %% Calculate mse
            fprintf('\tCalculate mmse...\n\n');
            y_hat = ((hr'*btheta*G*w*s)*LoS*LoS_PL+nLoS*LoS_PL+...
                (hd'*w*s)*Rayleigh_fading*nLoS_PL)+noise;
            s_hat = c*y_hat;
            if(t==1)
                mse = mean(power(abs(s_hat-s),2));
            elseif((mean(power(abs(s_hat-s),2))-mse)>epsilon)
                mse = mean(power(abs(s_hat-s),2));
            elseif((mean(power(abs(s_hat-s),2))-mse)<=epsilon)
                mse = mean(power(abs(s_hat-s),2));
                break;
            end
                       
            t = t+1;
        end
%% mode 2:The non-robust scheme
    elseif (mode==2)                            
        while true                              %Optimizes the system as if hr_hat, hd_hat and G_hat are perfect
            fprintf('The %d generation...\n',t);
            %% Optimization of c
            fprintf('\tOptimization of c...\n');
            A = (G_hat'*btheta'*hr_hat)*hr_hat'*btheta*G_hat+...          %According to the formula              
                hd_hat*hr_hat'*btheta*G_hat+G_hat'*btheta'*hr_hat*hd_hat'+...
                hd_hat*hd_hat'+variance*power(norm(hr_hat),2)*eye(M)+(variance*G_hat')*G_hat+...
                (N*variance*variance+variance)*eye(M);

            alpha = (G_hat'*btheta'*hr_hat+hd_hat);
        
            c = (w'*alpha)/(w'*A*w+noise_var);
        
            %% Optimization of w
            fprintf('\tOptimization of w...\n');
            w = (power(abs(c),2)*A)\(alpha*conj(c));  %Assume lambda(Lagrange mult.) = 0
            power(norm(w),2);
            if(power(norm(w),2)>P0)
                [lambda] = search(c,alpha,A,P0,M);
                w = (power(abs(c),2)*A+lambda*eye(M))\(alpha*conj(c));
            else
                lambda = 0;
                w = (power(abs(c),2)*A+lambda*eye(M))\(alpha*conj(c));
            end
            
            %% Optimization of theta
            fprintf('\tOptimization of btheta...\n');
            Phi = diag(hr_hat')*G_hat*w*c;
            d = hd_hat'*w*c;
            Q = Phi*Phi';
            q = Phi*(1-conj(d)');
            [new_v] = opt_v(Q,q,v);
            btheta = diag(new_v);
        
            %% Calculate mse
            fprintf('\tCalculate mmse...\n\n');
            y_hat = ((hr_hat'*btheta*G_hat*w*s)*LoS*LoS_PL+nLoS*LoS_PL+(hd_hat'*w*s)*Rayleigh_fading*nLoS_PL)+noise;
            s_hat = c*y_hat;
            if(t==1)
                mse = mean(power(abs(s_hat-s),2));
            elseif((mean(power(abs(s_hat-s),2))-mse)>epsilon)
                mse = mean(power(abs(s_hat-s),2));
                break;
            elseif((mean(power(abs(s_hat-s),2))-mse)<=epsilon)
                mse = mean(power(abs(s_hat-s),2));
                break;
            end            
            t = t+1;
        end
%% mode 3:The discrete phase shifts scheme
    elseif (mode==3)        
        while true                              
            fprintf('The %d generation...\n',t);
            %% Optimization of c
            fprintf('\tOptimization of c...\n');
            A = (G_hat'*btheta'*hr_hat)*hr_hat'*btheta*G_hat+...          %According to the formula              
                hd_hat*hr_hat'*btheta*G_hat+G_hat'*btheta'*hr_hat*hd_hat'+...
                hd_hat*hd_hat'+variance*power(norm(hr_hat),2)*eye(M)+(variance*G_hat')*G_hat+...
                (N*variance*variance+variance)*eye(M);

            alpha = (G_hat'*btheta'*hr_hat+hd_hat);
        
            c = (w'*alpha)/(w'*A*w+noise_var);
        
            %% Optimization of w
            fprintf('\tOptimization of w...\n');
            w = (power(abs(c),2)*A)\(alpha*conj(c));  %Assume lambda(Lagrange mult.) = 0
            power(norm(w),2);
            if(power(norm(w),2)>P0)
                [lambda] = search(c,alpha,A,P0,M);
                w = (power(abs(c),2)*A+lambda*eye(M))\(alpha*conj(c));
            else
                lambda = 0;
                w = (power(abs(c),2)*A+lambda*eye(M))\(alpha*conj(c));
            end
            %% Optimization of theta
            fprintf('\tOptimization of btheta...\n');
            Phi = diag(hr_hat')*G_hat*w*c;
            d = hd_hat'*w*c;
            Q = Phi*Phi';
            q = Phi*(1-conj(d)');
            [new_v] = opt_v(Q,q,v);
            [final_v] = discetet_phase(bit_of_phase,new_v);
            v = final_v;
            btheta = diag(v);
        
            %% Calculate mse
            fprintf('\tCalculate mmse...\n\n');
            y_hat = ((hr'*btheta*G*w*s)*LoS*LoS_PL+nLoS*LoS_PL+(hd'*w*s)*Rayleigh_fading*nLoS_PL)+noise;
            s_hat = c*y_hat;
            if(t==1)
                mse = mean(power(abs(s_hat-s),2));
            elseif((mean(power(abs(s_hat-s),2))-mse)>epsilon)
                mse = mean(power(abs(s_hat-s),2));
                break;
            elseif((mean(power(abs(s_hat-s),2))-mse)<=epsilon)
                mse = mean(power(abs(s_hat-s),2));
                break;
            end            
            t = t+1;
        end
%% mode 4:The scheme when IRS is not deployed(btheta==0)        
    else                    
        while true
            btheta = 0;
            fprintf('The %d generation...\n',t);
            %% Optimization of c
            fprintf('\tOptimization of c...\n');
            A = (G_hat'*btheta'*hr_hat)*hr_hat'*btheta*G_hat+...          %According to the formula              
                hd_hat*hr_hat'*btheta*G_hat+G_hat'*btheta'*hr_hat*hd_hat'+...
                hd_hat*hd_hat'+variance*power(norm(hr_hat),2)*eye(M)+(variance*G_hat')*G_hat+...
                (N*variance*variance+variance)*eye(M);

            alpha = (G_hat'*btheta'*hr_hat+hd_hat);
        
            c = (w'*alpha)/(w'*A*w+noise_var);
        
            %% Optimization of w
            fprintf('\tOptimization of w...\n');
            w = (power(abs(c),2)*A)\(alpha*conj(c));  %Assume lambda(Lagrange mult.) = 0
            power(norm(w),2);
            if(power(norm(w),2)>P0)
                [lambda] = search(c,alpha,A,P0,M);
                w = (power(abs(c),2)*A+lambda*eye(M))\(alpha*conj(c));
            else
                lambda = 0;
                w = (power(abs(c),2)*A+lambda*eye(M))\(alpha*conj(c));
            end
            %% Calculate mse
            fprintf('\tCalculate mmse...\n\n');
            y_hat = ((hr'*btheta*G*w*s)*LoS*LoS_PL+nLoS*LoS_PL+(hd'*w*s)*Rayleigh_fading*nLoS_PL)+noise;
            s_hat = c*y_hat;
            if(t==1)
                mse = mean(power(abs(s_hat-s),2));
            elseif((mean(power(abs(s_hat-s),2))-mse)>epsilon)
                mse = mean(power(abs(s_hat-s),2));
                break;
            elseif((mean(power(abs(s_hat-s),2))-mse)<=epsilon)
                mse = mean(power(abs(s_hat-s),2));
                break;
            end            
            t = t+1;
        end
    end
end
