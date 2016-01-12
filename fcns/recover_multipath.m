function [x_hat_group,x_hat_filt,x_hat,ite_num] ...
    = recover_multipath(...
    y,K,b,x_init,tau,del,ite,thres_eps)

% solve PCML by SPISTA
[x_hat,~,ite_num] = spista(...
    y,K,b,x_init,...
    tau,del,ite);
% filter residuals
x_hat_filt = x_hat;
x_hat_filt(x_hat_filt<max(x_hat_filt)/thres_eps) = 0;
% group depths
[ind_hat_group,val_hat_group]...
    = depth_grouping(x_hat_filt);
x_hat_group = zeros(size(x_hat_filt));
x_hat_group(round(ind_hat_group)) = val_hat_group;

end

