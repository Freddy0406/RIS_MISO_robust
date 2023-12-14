function [newlambda] = search(c,alpha,A,P0,M)
lambda_max = 100;
lambda_min = 0;
lambda = lambda_max/2;
i = 0;

while(1)
    w = (power(abs(c),2)*A+lambda*eye(M))\(alpha*conj(c));
    if(round( power( norm(w),2 ) , 5)<P0)                 %lambda
        lambda_max = lambda;
        lambda = lambda_min+(lambda_max-lambda_min)/2;
    elseif (round( power( norm(w),2 ) , 5)>P0)            %lambda
        lambda_min = lambda;
        lambda = lambda_min+(lambda_max-lambda_min)/2;
    elseif(round( power( norm(w),2 ) , 5) == P0)
        break;
    end
    i=i+1
    % if(i>1000)
    %     break;
    % end
end
    newlambda = lambda;
end