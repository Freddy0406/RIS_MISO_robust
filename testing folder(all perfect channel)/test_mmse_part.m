clc;
clear;
close all;

load('workspace');
t=1;

%% Start optimization    

while true
    %% Optimization of c
    hr_hat_norm = 0;
    for i = 1:length(hr_hat)
        hr_hat_norm = hr_hat_norm+abs(hr_hat(i));
    end
    hr_hat_norm = power(hr_hat_norm,2);

    A = G_hat'*btheta'*hr_hat*hr_hat'*btheta*G_hat+...
        hd_hat*hr_hat'*btheta*G_hat+...
        hd_hat*hd_hat'+variance*hr_hat_norm*eye(M)+variance*G_hat'*G_hat+...
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
    s_hat = c*(hr_hat'*btheta*G_hat+hd_hat')*w*s+noise;
    if(t==1)
        mse = mean(power(s_hat-s,2));
    else
        break;
    end
   
    t = t+1;
end



