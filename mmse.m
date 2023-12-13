function [mse,t]=mmse(N,variance,P0,mode)
%% Initialization parameter
    %N Reflecting element amount per RIS
    %P0 Maximum transmit power of BS
    %variance var_g=var_r=var_d=var
    M = 4;                          %BS antenna amount
    s = randn(1,1,"like",1i);       %Symbol for user(zero-meam, unit power=>variance=1)
    L0 = -30;                       %Path loss dB           
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


    % Initialize w
    w = zeros(M,1);

    for i = 1:M
        w(i) = sqrt(P0)/sqrt(M);
    end

    % Initialize theta
    min = 0;
    max = 2*pi;
    v = exp(1i*(min+rand(N,1)*(max-min)));
    btheta = diag(v);

    % delta_G, delta_hd, delta_hr (AP-IRS, IRS-User, AP-User link)
    
    delta_G = sqrt(variance)*randn(N,M,"like",1i);
    delta_hd = sqrt(variance)*randn(M,1,"like",1i);
    delta_hr = sqrt(variance)*randn(N,1,"like",1i);
    
    % Emulate estimation channel

    G = vertcat(eye(M),zeros((N-M),M));      %Ideal channel
    hd = ones(M,1);                         
    hr = ones(N,1);
    
    G_hat = G-delta_G;
    hd_hat = hd-delta_hd;
    hr_hat = hr-delta_hr;

    t=1;
    %% Start optimization    
    if (mode==1)            %%mode 1:The proposed robust design
        while true
            fprintf('The %d generation...\n',t);
            %% Optimization of c
            fprintf('\tOptimization of c...\n');
            hr_hat_norm = 0;
            for i = 1:length(hr_hat)
                hr_hat_norm = hr_hat_norm+abs(hr_hat(i));
            end
            hr_hat_norm = power(hr_hat_norm,2);
        
            A = (G_hat'*btheta'*hr_hat)*hr_hat'*btheta*G_hat+...          %According to the formula              
                hd_hat*hr_hat'*btheta*G_hat+G_hat'*btheta'*hr_hat*hd_hat'+...
                hd_hat*hd_hat'+variance*hr_hat_norm*eye(M)+(variance*G_hat')*G_hat+...
                (N*variance*variance+variance)*eye(M);
            alpha = (G_hat'*btheta'*hr_hat+hd_hat);
        
            c = (w'*alpha)/(w'*A*w+noise_var);          %Wiener Filter
        
            %% Optimization of w
            fprintf('\tOptimization of w...\n');
            w = (power(abs(c),2)*A)\(alpha*conj(c));  %Assume lambda(Lagrange mult.) = 0
            power(norm(w),2);
            if(power(norm(w),2)>P0)
                [lambda] = search(c,alpha,A,P0,M);
                w = power((power(abs(c),2)*A+lambda*eye(M)),-1)*(alpha*conj(c));
            else
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
            y_hat = ((hr'*btheta*G*w*s)*LoS_PLE+(hd'*w*s)*nLoS_PLE)+noise;
            s_hat = c*y_hat;
            if(t==1)
                mse = mean(power(abs(s_hat-s),2));
            elseif((mean(power(abs(s_hat-s),2))-mse)>epsilon)
                mse = mean(power(abs(s_hat-s),2));
            elseif((mean(power(abs(s_hat-s),2))-mse)<=epsilon)
                break;
            end
                       
            t = t+1;
        end

    elseif (mode==2)                            %mode 2:The non-robust scheme
        while true
            %% Optimization of c
            hr_hat_norm = 0;
            for i = 1:length(hr_hat)
                hr_hat_norm = hr_hat_norm+abs(hr_hat(i));
            end
            hr_hat_norm = power(hr_hat_norm,2);
        
            A = (G_hat'*btheta'*hr_hat)*hr_hat'*btheta*G_hat+...
                hd_hat*hr_hat'*btheta*G_hat+...
                hd_hat*hd_hat'+variance*hr_hat_norm*eye(M)+(variance*G_hat')*G_hat+...
                (N*variance*variance+variance)*eye(M);
            alpha = (G_hat'*btheta'*hr_hat+hd_hat);
        
            c = (w'*alpha)/(w'*A*w+noise_var);
        
            %% Optimization of w
            w = power((power(abs(c),2)*A),-1)*(alpha*conj(c));
        
            if(power(norm(w,1),2)>P0)
                [lambda] = search(c,alpha,A,P0,M);
                w = power(((power(abs(c),2)*A+lambda*eye(M))),-1)*(alpha*conj(c));
            else
            end
            
            %% Optimization of theta
            Phi = diag(hr_hat')*G_hat*w*c;
            d = hd_hat'*w*c;
            Q = Phi*Phi';
            q = Phi*(1-conj(d)');
            [new_v] = opt_v(Q,q,v);
            btheta = diag(new_v);
        
            %% Calculate mse
            y_hat = ((hr_hat'*btheta*G_hat*w*s)*LoS_PLE+(hd_hat'*w*s)*nLoS_PLE)+noise;
            s_hat = c*y_hat;
            if(t==1)
                mse = mean(power(abs(s_hat-s),2));
            elseif((mean(power(abs(s_hat-s),2))-mse)>power(10,-4))
                mse = mean(power(abs(s_hat-s),2));
                break;
            elseif((mean(power(abs(s_hat-s),2))-mse)<=power(10,-4))
                break;
            end
           
            t = t+1;
        end

    elseif (mode==3)        %mode 3:The discrete phase shifts scheme
    else                    %mode 4:The scheme when IRS is not deployed(btheta==0)

    end
end
