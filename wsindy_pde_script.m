%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WSINDy_PDE: script for recoverying PDE systems
%%%%%%%%%%%% 
%%%%%%%%%%%% pde_num selects a PDE system from the list pde_names
%%%%%%%%%%%% noise_ratio sets the signal-to-noise ratio (L2 sense)
%%%%%%%%%%%% 
%%%%%%%%%%%% Copyright 2020, All Rights Reserved
%%%%%%%%%%%% Code by Daniel A. Messenger
%%%%%%%%%%%% For Paper, "Weak SINDy for Partial Differential Equations"
%%%%%%%%%%%% by D. A. Messenger and D. M. Bortz

%% Load data

clc; 
% clear all; 
% close all;

pde_num = 3;
pde_names = {'burgers.mat','KdV.mat','KS.mat','NLS.mat','Sine_Gordon.mat','rxn_diff.mat','Nav_Stokes.mat','porous.mat','bw'};
load(['datasets/',pde_names{pde_num}])
% load(['/home/danielmessenger/Dropbox/Boulder/research/data/WSINDy_PDE/datasets/',pde_names{pde_num}])

%%

U_obs = U_exact;
xs_obs = xs;
dims = size(U_obs{1});
dim = length(dims);
n = length(U_obs);

%% Subsample data (if desired)

coarsen_data = [ones(dim,2) dims'];
coarsen_data(:,2) = 1;
%%% set row d to [i inc f] to subsample dth coordinate to start at index i,
%%% end at index f, and skip every inc gridpoint. e.g:
%%% coarsen_data(:,2) = 3; coarsens spatiotemporal grid by factor of 3 in each coordinate

[xs_obs,U_obs] = subsamp(xs_obs,U_obs,coarsen_data,dims);
dims = cellfun(@(x) length(x), xs_obs);

%% Add noise

sigma_NR = 0;
noise_dist = 0; 
noise_alg = 0;
rng('shuffle');
rng_seed = rng().Seed;
 
rng(rng_seed);
[U_obs,noise,snr,sigma] = gen_noise(U_obs,sigma_NR,noise_dist,noise_alg,rng_seed,0);

%% Set hyperparameters 

%---------------- weak discretization

phi_class = 1;          % piecewise poly = 1, Gaussian = 2

s_x = 5;  %max(floor(length(xs_obs{1})/50),1);
s_t = 5;  %max(floor(length(xs_obs{end})/50),1);

m_x = 10;               % test function in x, 1/2 support in grid points 
m_t = 10;               % test function in t, 1/2 support in grid points 
p_x = 4;               % test function in x, poly degree
p_t = 2;               % test function in t, poly degree

tauhat = 4;             % if 0, set m_x,m_t,p_x,p_t manually; where tau = if 1, use cornerpoint
tau = 10^-10;

toggle_scale = 2;

%---------------- model library

max_dx = 5;
max_dt = 1;
polys = 0:5;
trigs = [];
use_all_dt = 0;
use_cross_dx = 0;

custom_add = [];
custom_remove = [];
% lhs = [];

%---------------- find test function hyperparams using Fourier spectrum of U

if tauhat > 0
    tauhat_inds = 1:n;
    [m_x,m_t,p_x,p_t,sig_est,corners_all] = findcorners(cellfun(@(x) x.^1, U_obs(tauhat_inds), 'uni',0),xs_obs,tau,tauhat,max_dx,max_dt,phi_class);
else
    m_x = min(m_x,floor((length(xs_obs{1})-1)/2));
    m_t = min(m_t,floor((length(xs_obs{end})-1)/2));
end
tols = [-p_x -p_t];

%% Build Library

[axi,tags_pde,lib_list,pdx_list,lhs_ind,Cfs_x,Cfs_t,dx,dt,p_x,p_t,sub_inds,scales,M_full,Theta_pdx] = wsindy_pde_fun(U_obs,xs_obs,true_nz_weights,...
    lhs,max_dx,max_dt,polys,trigs,custom_add,custom_remove,use_all_dt,use_cross_dx,...
    toggle_scale,m_x,m_t,s_x,s_t,tols,phi_class);

%% Solve Sparse Regression Problem

%---------------- regularized least squares solve
lambda = 10.^(linspace(-4,0,100));
% lambda = [];
gamma = 0;
maxits = Inf;
sparsity_scale = 1;                     % 0, enforce sparsity on original data; 1, enforce on rescaled data

[W,G,b,resid,dW,its_all,thrs_EL,M,lambda_hat,lossvals,ET_wsindy,tags_pde_G,lib_list_G] = wsindy_pde_solve(lambda,gamma,Theta_pdx,lhs_ind,axi,M_full,maxits,tags_pde,lib_list,sparsity_scale);

%%

print_loc = 1;

toggle_plot_basis_fcn = 1;
toggle_plot_sol =  1;
toggle_plot_loss = 1;
toggle_plot_fft = 1;

get_results;
display_results;

% %% Display results
% 
% print_loc = 1;
% toggle_plot_basis_fcn = 0; 
% toggle_plot_sol = 0; 
% toggle_plot_loss = 1; 
% 
% if print_loc~=0
%     [str_wsindy,Tpr]  = print_results(W,G,resid,dW,print_loc,dims,polys,trigs,max_dx,max_dt,lambda_hat,gamma,lhs_ind,tags_pde,m_x,m_t,p_x,p_t,sub_inds,scales,ET_wsindy,its_all,sigma_NR,sigma,axi);
% end
% 
% if toggle_plot_basis_fcn
%     figure(1);
%     plot_basis_fcn(Cfs_x,Cfs_t,m_x,dx,m_t,dt,max_dx,max_dt,pdx_list,[1000 5 1400 700],scales(end-n:end));
% end
% 
% if toggle_plot_sol>0
%     U = U_obs{1};
%     figure(2); %set(gcf,'units','normalized','outerposition',[0.7 0.5 0.3 0.45])
%     set(gca, 'TickLabelInterpreter','latex','fontsize',14)
%     xlabel('$x$','interpreter','latex','fontsize',14)
%     ylabel('$t$','interpreter','latex','fontsize',14)
%     if dim==2
%         surf(1:length(xs_obs{1}),1:length(xs_obs{2}),U', 'EdgeColor','none')
%         colorbar
%         view([0 90])       
%         zlabel('$u$','interpreter','latex','fontsize',14)
%         set(gca, 'TickLabelInterpreter','latex','fontsize',14)
%     elseif dim==3
%         for j=1:toggle_plot_sol:length(xs_obs{3})
%             surf(xs_obs{1},xs_obs{2},U(:,:,j)', 'EdgeColor','none')
%             zlim([0 max(U(:))])
%             colormap(turbo)
%             view([0 90])
%             caxis([min(U(:)) max(U(:))])
%             drawnow
%         end
%     end
% end
% 
% if ~isempty(lossvals) && toggle_plot_loss
%     figure(3); %set(gcf,'units','normalized','outerposition',[0.7 0.05 0.3 0.45])
%     loglog(lossvals(2,:),lossvals(1,:),'o-')
%     xlabel('$\lambda$','interpreter','latex','fontsize',14)
%     ylabel('$\mathcal{L}$','interpreter','latex','fontsize',14)
%     set(gca, 'TickLabelInterpreter','latex','fontsize',14)
%     xtick =10.^linspace(log10(lossvals(2,1)),log10(lossvals(2,end)),4);
%     xticks(xtick)
%     xticklb = num2str(log10(xtick)');
%     xticklabels(strcat('$10^{',xticklb(:,1:min(4,end)),'}$'))
% end
% 
