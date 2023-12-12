function [new_v] = opt_v(Q,q,v)

iteration = 0;

while true

    [V,D] = eig(Q);
    lambda_Q = max(D,[],"all");
    
    u = (Q-lambda_Q*eye(length(Q)))*v-q;
    
    next_v = -exp(1i*angle(u));

    if(sum(next_v-v)==0)
        break;
    elseif(iteration>1000)
        break;
    else
    end
    iteration = iteration+1;
    v = next_v;
end   
    new_v = next_v;
end