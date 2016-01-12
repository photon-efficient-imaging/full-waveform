function [ x_hat, x_hat_history, c ] = spista(...
    y,K,b,x_init,...
    tau,del,max_ite)

is_history = 0;
x_hat_history = [];

x_hat = x_init;
x_hat_prev = zeros(size(x_hat));
c = 0;
errs = 1e2;
while(errs>del)
    nabla_x = K(1 - (y./(K(x_hat)+b+eps)));
    %nabla_x = K(K(x_hat)-y);

    x_hat = x_hat - nabla_x;
    % 2. soft threshold
    x_hat = max(x_hat - tau, 0.0);
    errs = norm(x_hat-x_hat_prev);
    x_hat_prev = x_hat;
    c = c+1;
    if(c==max_ite)
        break;
    end
    
    if(is_history)
       x_hat_history = [x_hat_history, x_hat]; 
    end
end
end

