function [mse,t]=mmse(mode)

    load('workspace');
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
        
            A = G_hat'*btheta'*hr_hat*hr_hat'*btheta*G_hat+...          %According to the formula              
                hd_hat*hr_hat'*btheta*G_hat+G_hat'*btheta'*hr_hat*hd_hat'+...
                hd_hat*hd_hat'+variance*hr_hat_norm*eye(M)+variance*G_hat'*G_hat+...
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
            elseif((mean(power(s_hat-s,2))-mse)>power(10,-4))
                mse = mean(power(s_hat-s,2));
                break;
            elseif((mean(power(s_hat-s,2))-mse)<=power(10,-4))
                break;
            end
           
            t = t+1;
        end
    else                    %mode 4:The scheme when IRS is not deployed(btheta==0)

    end
end
