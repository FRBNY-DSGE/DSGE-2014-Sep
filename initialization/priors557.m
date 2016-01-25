% Define Prior parameters
% pshape is 1: BETA(mean,stdd)
%           2: GAMMA(mean,stdd)
%           3: NORMAL(mean,stdd)
%           4: INVGAMMA(s^2,nu)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% loose lambda_f prior %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function prior = priors557()

prior = zeros(100,3);
nantpad = 20;



prior(1,:) = [.33   .02  1];     %% alp - beta
prior(2,:) = [.75    .1  1];  %% zeta_p - beta
prior(3,:) = [.5    .28 1];  %% iota_p - betaa
prior(4,:) = [.025  .01 1]; %% del - beta
prior(5,:) = [963.55 10 3];   %% Ztil_0 - normal 963.55
prior(6,:) = [.5    .25 2];  %% Bigphi - gamma
prior(7,:) = [4     1.5 2];  %% s2 - normal
prior(8,:) = [.7    .05 1];  %% h - beta
prior(9,:) = [.2    .10 2];  %% a2 - gamma
prior(10,:) = [2    .75 2];  %% nu_l - gamma
prior(11,:) = [2    .75 2];  %% nu_m - gamma
prior(12,:) = [.75   .1  1];  %% zeta_w - beta
prior(13,:) = [.5   .28  1]; %% iota_w - beta
prior(14,:) = [.3   .5  2];  %% law - gamma
prior(15,:) = [1.5   1  2];  %% rstar - gamma
prior(16,:) = [2  .25  2];  %% psi1 - gamma
prior(17,:) = [.2   .1  2];  %% psi2 - gamma
prior(18,:) = [.5   .2  1];  %% rho_r - beta
prior(19,:) = [2  0.25  3];  %% pistar - normal
prior(20,:) = [.03 .01  1];  %% Fom - beta
prior(21,:) = [2    .5   2];  %% st st spread - gamma
prior(22,:) = [.05 .02  1];  %% zeta_sp - beta
prior(23,:) = [.99 .002 1];  %% gammstar - beta
prior(24,:) = [0 1  3];  %% NEW_5 - normal

npara = 24; 

%% exogenous processes - level
prior(npara+1,:) = [2.75 .5 2];  %% gam - gamma
prior(npara+2,:) = [0     5 3];     %% Wadj - normal
prior(npara+3,:) = [.1   .1  2];    %% chi - gamma
prior(npara+4,:) = [.15  .1  2];    %% laf - gamma
prior(npara+5,:) = [.3   .1  2]; %% gstar - beta
prior(npara+6,:) = [253.5  5 3];     %% Ladj - normal

npara = npara+6;    

%% exogenous processes - autocorrelation
prior(npara+1,:) = [.4   .25  1];   %% rho_z - beta
prior(npara+2,:) = [.75  .15  1];   %% rho_phi - beta
prior(npara+3,:) = [.75  .15  1];   %% rho_chi - beta
prior(npara+4,:) = [.75  .15  1];   %% rho_laf - beta
prior(npara+5,:) = [.75  .15  1];   %% rho_mu - beta
prior(npara+6,:) = [.75  .15  1];   %% rho_b - beta
prior(npara+7,:) = [.75  .15  1];   %% rho_g - beta
prior(npara+8,:) = [.75  .15  1];   %% rho_sigw - beta
prior(npara+9,:) = [.75  .15  1];   %% rho_mue - beta
prior(npara+10,:) = [.75  .15  1];   %% rho_gamm - beta

npara = npara+10;    

%% exogenous processes - standard deviation
prior(npara+1,:) = [.3   4.00    4]; %% sig_z;
prior(npara+2,:) = [3    4.00    4]; %% sig_phi;
prior(npara+3,:) = [.75  4.00    4]; %% sig_chi;
prior(npara+4,:) = [.2  4.00    4]; %% sig_laf;
prior(npara+5,:) = [.75  4.00    4]; %% sig_mu;
prior(npara+6,:) = [.75  4.00    4]; %% sig_b;
prior(npara+7,:) = [.5  4.00    4]; %% sig_g;
prior(npara+8,:) = [.20  4.00    4]; %% sig_r;
prior(npara+9,:) = [.20/4  4.00    4]; %% sig_sigw;
prior(npara+10,:) = [.20/4  4.00    4]; %% sig_mue;
prior(npara+11,:) = [.01  4.00    4]; %% sig_gamm;

npara = npara+11;

for i = 1:nantpad
    eval(strcat('prior(npara +',num2str(i),',:) = [.2  4.00    4];'));
end

npara = npara+nantpad;

prior = prior(1:npara,:);
