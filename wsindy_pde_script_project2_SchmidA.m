% APPM 5720-038, Spring 2022
% Project 2
% Abby Schmid

% script to recover PDE systems
% modified from wsindy_pde_script from the paper
% "Weak SINDy for Partial Differential Equations"
% by D. A. Messenger and D. M. Bortz

%% Choose dataset and data mode

clc; 

% choose data set, must be 1 or 2
% 1 = 01-25, sine wave input 
% 2 = 02-09, 5 cycle burst input  
data_set = 1;

% choose data mode, must be 'velocity' or 'displacement'
% velocity = experimentally measured velocity values
% displacement = numerically integrated experimental data
data_mode = 'displacement'; 

%% Load data

if data_set == 1
    velocity_data = importdata('data\velocity-clean-Al1-01-25.mat');
elseif data_set == 2
    velocity_data = importdata('data\velocity-clean-Al1-02-09.mat');
else
    error('data_set must be 1 or 2')
end 

[x_pts, t_pts] = size(velocity_data);

% integer steps
%x_vec = 0:x_pts-1;
%t_vec = 0:t_pts-1;

% physically relevant steps
x_vec = linspace(0, .1, x_pts); %beam length=.1m
t_vec = (1e-6)*(0:t_pts-1); %frequency=1e6, 1e-6 s

displacement_data = integrate_data(velocity_data, t_vec);

% convert to cell data structure
if strcmp(data_mode, 'velocity')
    U_obs = num2cell(velocity_data,[1 2]);
elseif strcmp(data_mode, 'displacement')
    U_obs = num2cell(displacement_data,[1 2]);
else
    error('data_mode must be "velocity" or "displacement"')
end

dims = size(U_obs{1});
xs_obs = num2cell(x_vec, [2 1]);
xs_obs{2} = num2cell(t_vec, 1);
xs_obs{2} = cell2mat(xs_obs{2});
dim = length(dims);
n = length(U_obs);

%% Add noise
% keep at 0, working with experimental data 
sigma_NR = 0.0;
noise_dist = 0; 
noise_alg = 0;
rng('shuffle');
rng_seed = rng().Seed;
 
rng(rng_seed);
[U_obs,noise,snr,sigma] = gen_noise(U_obs,sigma_NR,noise_dist,noise_alg,rng_seed,0);

%% set plots

toggle_plot_basis_fcn = 1;
toggle_plot_sol =  1;
plotgap = 3;
toggle_plot_loss = 1;
toggle_plot_fft = 1;

%% Set hyperparameters 

%---------------- weak discretization
%%% phi_class = 1 for piecewise polynomial test function, 2 for Gaussian
phi_class = 1;          

%%% set convolution query point spacing:
%s_x = 5;  
%s_t = 5; 
s_x = max(floor(length(xs_obs{1})/50),1);
s_t = max(floor(length(xs_obs{end})/50),1);

%%% set reference test function parameters using spectrum of data:
%%% if tauhat<=0, explicit vals for m_x,m_t,p_x,p_t used. 
% tauhats chosen to have the red line (FT of the test function) 
% go through the center of the yellow dot (corner of the data FT spectrum) 
% for the given data set and phi_class 
if data_set == 1
    if phi_class == 1
        tauhat = 2;
    else %phi_class==2
        tauhat = 2.25;
    end
else %data_set == 2
    tauhat = 1.25;
end

tau = 10^-10;

%%% set reference test function parameters explicitly:
m_x = 15;
m_t = 15;
p_x = 11;
p_t = 7;

%%% toggle rescale state variables and spatiotemporal coordinates
toggle_scale = 2;

%---------------- model library
if strcmp(data_mode, 'velocity')
    max_dt = 1; %up to acceleration 
    % lhs = time derivative of velocity 
    lhs = [1 0 1];
else %strcmp(data_mode, 'displacement')
    max_dt = 2; %up to acceleration 
    % lhs = second time derivative of displacement 
    lhs = [1 0 2];
end

max_dx = 4;
polys = 0:4;
trigs = [];
use_all_dt = 1;
use_cross_dx = 0;
custom_add = [];
custom_remove = {}; %{@(mat,lhs) find(all([mat(:,3)==0 ~ismember(mat,lhs,'rows')],2))};

% set lhs
% 1 x (n+dim) row vector, [p1 ... pn ... d1 ... ddim], n=num state variables, dim = dimensions
% if strcmp(data_mode, 'velocity')
%     % lhs = time derivative of velocity 
%     lhs = [1 0 1];
% else %strcmp(data_mode, 'displacement')
%     % lhs = second time derivative of displacement 
%     lhs = [1 0 2];
% end
 
true_nz_weights = {};

%% Build Linear System

%---------------- find test function hyperparams using Fourier spectrum of U
if tauhat > 0
    tauhat_inds = 1:n;
    [m_x,m_t,p_x,p_t,sig_est,corners] = findcorners(cellfun(@(x) x.^1, U_obs(tauhat_inds), 'uni',0),xs_obs,tau,tauhat,max_dx,max_dt,phi_class);
else
    m_x = min(m_x,floor((length(xs_obs{1})-1)/2));
    m_t = min(m_t,floor((length(xs_obs{end})-1)/2));
end
tols = [-p_x -p_t];

%---------------- build linear system
[axi,tags_pde,lib_list,pdx_list,lhs_ind,Cfs_x,Cfs_t,dx,dt,p_x,p_t,sub_inds,scales,M_full,Theta_pdx] = wsindy_pde_fun(U_obs,xs_obs,true_nz_weights,...
    lhs,max_dx,max_dt,polys,trigs,custom_add,custom_remove,use_all_dt,use_cross_dx,...
    toggle_scale,m_x,m_t,s_x,s_t,tols,phi_class);

%% Solve Sparse Regression Problem

lambda = 10.^(linspace(-4,0,50));
gamma = 0;
maxits = Inf;

%%% sparsity_scale =  0 enforces sparsity on original data; = 1 enforces on rescaled data
sparsity_scale = 0;                     

[W,G,b,resid,dW,its_all,thrs_EL,M,lambda_hat,lossvals,ET_wsindy,tags_pde_G,lib_list_G] = wsindy_pde_solve(lambda,gamma,Theta_pdx,lhs_ind,axi,M_full,maxits,tags_pde,lib_list,sparsity_scale);

%% Display results

print_loc = 1;
get_results;
display_results;
