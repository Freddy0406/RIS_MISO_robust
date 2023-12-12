function [newlambda] = bisection(c,alpha,A,P0,M)
iteration = 0;
lambda_max = 10;
lambda_min = 0;
lambda = lambda_max;


while(1)
    w = power((power(abs(c),2)*A+lambda*eye(M)),-1)*(alpha*conj(c));

    if(power(norm(w,1),2)<P0)
        lambda_max = lambda;
    elseif (power(norm(w,1),2)>P0)
        lambda_min = lambda;
    elseif(power(norm(w,1),2) == P0)
        break;
    end
    lambda = (lambda_max-lambda_min)/2;

    iteration=iteration+1;
    if(iteration>1000)
        break;
    end
end
    newlambda = lambda;
end