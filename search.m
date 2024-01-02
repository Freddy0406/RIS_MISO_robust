function [newlambda] = search(c,alpha,A,P0,M)
lambda_max = 10000;
lambda_min = 0;
lambda = lambda_max;
i = 0;

while(1)
    w = (power(abs(c),2)*A+lambda*eye(M))^-1*(alpha*conj(c));
    if(round( power( norm(w),2 ) , 5)<P0)                 %lambda needs smaller
        lambda_max = lambda;
        lambda = lambda_min+(lambda_max-lambda_min)/2;
    elseif (round( power( norm(w),2 ) , 5)>P0)            %lambda needs larger
        lambda_min = lambda;
        lambda = lambda_min+(lambda_max-lambda_min)/2;
    elseif(round( power( norm(w),2 ) , 5) == P0)
        break;
%     elseif (i>100000)
%         break;
    end
    i=i+1;
end
    newlambda = lambda;
end