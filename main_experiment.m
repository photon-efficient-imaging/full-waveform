% Computational single-photon multi-depth imaging
% D.Shin, F.Xu, F.N.C.Wong, J.H.Shapiro, V.K.Goyal

clc; close all; clear;

% load data and set imaging parameters
load('two_layer');
[nr,nc] = size(photon_times);
kk = 91;
sigs = 35;
xxx = -kk:kk;
kernel = exp(-(xxx.^2)./(2*sigs^2));
sumkernel = sum(kernel);
K = @(inputs) conv(inputs,kernel,'same');
K_norm = @(inputs) conv(inputs,kernel,'same')/sumkernel;
v = 3000:7000;
% peak is 1 for narrow-width normal-approximated pulse
ff = @(x) exp(-abs(x).^2/(2*sigs^2)); 
nn = length(v);
t1 = 1:nn;
S = zeros(nn);
for i=1:nn    
    s = ff(t1-t1(i)); s = s/max(s); S(:,i) = s';
end
% relevant time-ranges for plotting
min_r_1 = 4200; max_r_1 = 4900;
min_r_2 = 5900; max_r_2 = 6500;
%
run_gmm = 0; % run gmm by setting this to 1
% In order to run code with run_gmm=1
% download EM-based MoG method (written by S. Chen) from
% http://www.mathworks.com/matlabcentral/fileexchange
%       /26184-em-algorithm-for-gaussian-mixture-model/content/emgm/emgm.m
% and place /emgm in local directory
addpath(genpath([pwd '/emgm']));

% some parameters
% [1] gmm
s_gmm = 2; % estimate two reflectors
% [2] ours
val_bgd = 0.02;
val_tau = 0.1;
val_del = 1e-4;
val_max_ite = 100;
val_thres_eps = 3;
addpath(genpath([pwd '/fcns']));
D_sparse = cell(nr,nc); D_gmm = cell(nr,nc);
D_sparse_1 = zeros(nr,nc); D_sparse_2 = zeros(nr,nc);
D_gmm_1 = zeros(nr,nc); D_gmm_2 = zeros(nr,nc);
fprintf('# Running reconstruction algorithms ...\n');
for i=1:nr
    fprintf(['row ' num2str(i) ' out of ' num2str(nr) '\n']);
    for j=1:nc
        dats_fin = photon_times{i,j};
        y = hist(dats_fin,v);
        y = y';
        % [1] gmm estimator
        if(run_gmm)
        ind_gmm = [];
        if(length(find(y))~=1)            
            samples = hist2sample(y);
            [vv,gmm_struct] = emgm(samples,s_gmm);
            ind_gmm = round(gmm_struct.mu);
            a_gmm = gmm_struct.weight;
            sol_gmm = v(ind_gmm);
            D_gmm{i,j} = sol_gmm;           
            
            m_1 = (sol_gmm>min_r_1)&(sol_gmm<max_r_1);
            m_2 = (sol_gmm>min_r_2)&(sol_gmm<max_r_2);
            if(~isempty(sol_gmm(m_1)))
                sss = a_gmm(m_1);
                sols = sol_gmm(m_1);
                [vals,inds] = max(sss);
                D_gmm_1(i,j) = sols(inds);
            end
            if(~isempty(sol_gmm(m_2)))
                sss = a_gmm(m_2);
                sols = sol_gmm(m_2);
                [vals,inds] = max(sss);
                D_gmm_2(i,j) = sols(inds);
            end
        end
        end
                
        % [2] sparsity pursuit        
        x_init = conv(y,kernel,'same')/sum(kernel);
        [x_hat_group,x_hat_filt,x_hat,ite_SPISTA] ...
            = recover_multipath(...
            y,K_norm,val_bgd,x_init,...
            val_tau,val_del,val_max_ite,val_thres_eps);        
        i_hat_group = find(x_hat_group);
        a_hat_group = x_hat_group(i_hat_group);
        sol_sparse = v(i_hat_group);
        D_sparse{i,j} = sol_sparse;
        m_1 = (sol_sparse>min_r_1)&(sol_sparse<max_r_1);
        m_2 = (sol_sparse>min_r_2)&(sol_sparse<max_r_2);
        if(~isempty(sol_sparse(m_1)))
            sss = a_hat_group(m_1);
            sols = sol_sparse(m_1);
            [vals,inds] = max(sss);
            D_sparse_1(i,j) = sols(inds);
        end
        if(~isempty(sol_sparse(m_2)))
            sss = a_hat_group(m_2);
            sols = sol_sparse(m_2);
            [vals,inds] = max(sss);
            D_sparse_2(i,j) = sols(inds);
        end

    end
end
fprintf('done! \n')

ranges = [5900,6180];
load('two_layer_obj')
figure;
subplot(221); imagesc(T_first); 
axis image; title({'calibrated depth','of scatterer'})
colormap('jet')
subplot(222); imagesc(T_second,ranges); 
axis image; title({'calibrated depth','of mannequin behind'})
D1 = D_gmm_2; D1(D1==0) = inf;
subplot(223); imagesc(D1,ranges); axis image;
title({'MoG estimate of','mannequin depth'});
D2 = D_sparse_2; D2(D2==0) = inf;
subplot(224); imagesc(D2,ranges); axis image;
title({'Proposed estimate of','mannequin depth'});

