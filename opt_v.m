function [new_v] = opt_v(Q,q,v)
i = 0;
while true

    D = eig(Q);                 %Eigenvalue
    lambda_Q = max(D,[],"all"); %Find maximum of Eigenvalue
    
    u = (Q-lambda_Q*eye(length(Q)))*v-q;
    
    next_v = -exp(1i*angle(u));

    if(sum(next_v-v)<power(10,-20))
        break;
    else
    end
    v = next_v;  
    i=i+1;
    if(i>40000)
        break;
    end
end   
    new_v = next_v;
end