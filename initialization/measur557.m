%% 
% measur557.m:
% solution to DSGE model - delivers transition equation for the 
% state variables S_t
% transition equation: S_t = TC+TTT S_{t-1} +RRR eps_t, where var(eps_t) = QQ
% define the measurement equation: X_t = ZZ S_t +D+u_t
% where u_t = eta_t+MM* eps_t with var(eta_t) = EE
% where var(u_t) = HH = EE+MM QQ MM', cov(eps_t,u_t) = VV = QQ*MM'

function [ZZ,DD,DDcointadd,QQ,EE,MM,retcode] = ...
  measur557(TTT,RRR,valid,para,nvar,nlags,mspec,npara,coint,cointadd,nant);

%Exit with return code of 0 if inputs not valid.
retcode = 1;
if valid < 1;
    retcode = 0;
    ZZ = [];
    DD = [];
    QQ = [];
    EE = [];
    MM = [];    
    DDcointadd = [];
    return
end

nstate = size(TTT,1);
DDcointadd = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% step 1: assign names to the parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% Parameters 
% [alp,zeta_p,iota_p,del,ups,Bigphi,s2,h,a2,nu_l,nu_m,zeta_w,iota_w,law,rstar,psi1,psi2,rho_r,pistar,...
%           Fom,sprd,zeta_spb,gammstar,NEW_5,...
%           gam,Lstar,chi,laf,gstar,Ladj,...
%           rho_z,rho_phi,rho_chi,rho_laf,rho_mu,rho_b,rho_g,rho_sigw,rho_mue,rho_gamm,...
%           sig_z,sig_phi,sig_chi,sig_laf,sig_mu,sig_b,sig_g,sig_r,sig_sigw,sig_mue,sig_gamm,...
%           bet,zstar,phi,istokbarst,rkstar,cstar,ystar,istar,kstar,kbarstar,Rstarn,wstar,wadj,mstar,...
%           zeta_nRk, zeta_nR, zeta_nqk, zeta_nn, zeta_nmue, zeta_spmue, zeta_nsigw, zeta_spsigw, ...
%           vstar, nstar,Rkstar, sig_r_ant] = getpara00_557(para);
getPara_script;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% step 2: assign names to the columns of GAM0, GAM1 -- state variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eval(strcat('states',num2str(mspec)));

%% additional states
y_t1 = n_end+n_exo+n_exp+1;
%c_t1 = n_end+n_exo+n_exp+2;
%i_t1 = n_end+n_exo+n_exp+3;
%w_t1 = n_end+n_exo+n_exp+4;
%m_t1 = n_end+n_exo+n_exp+5;

    if nstate ~= (n_end+n_exo+n_exp+1)

        retcode = 0;

        yyyyd = zeros(nvar,nvar);
        xxyyd = zeros(1+nlags*nvar,nvar);
        xxxxd = zeros(1+nlags*nvar,1+nlags*nvar);

        disp('\n\n number of states does not match in vaprio\n');
        return
    end
if ~exist('nant','var')
    nvar = 7;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% step 3: assign measurement equation : X_t = ZZ S_t +D+u_t
%% where u_t = eta_t+MM* eps_t with var(eta_t) = EE
%% where var(u_t) = HH = EE+MM QQ MM', cov(eps_t,u_t) = VV = QQ*MM'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create system matrices for state space model
ZZ = zeros(nvar+coint,nstate);
%% constant
DD = zeros(nvar+coint,1);
%% cov(eps_t,u_t) = VV
MM = zeros(nvar+coint,nex);   
%% var(eta_t) = EE
EE = zeros(nvar+coint,nvar+coint);
%% var(eps_t) = QQ
QQ =  zeros(nex,nex);

%% Output growth - Quarterly Annualized! 
ZZ(1,y_t) = 4;
ZZ(1,y_t1) = -4;
ZZ(1,z_t) = 4;
DD(1) = 400*(gam+(alp*log(ups)/(1-alp)));

%% Hours
ZZ(2,L_t) = 1;
DD(2) = 100*log(Ladj);%+log(Lstar);

%% Labor Share
ZZ(3,L_t) = 1;
ZZ(3,w_t) = 1;
ZZ(3,y_t) = -1;
DD(3) = 100*log((1-alp)/(1+laf));

%% Inflation
ZZ(4,pi_t) = 4;
DD(4) = 400*log(pistar);

%% Nominal interest rate
ZZ(5,R_t) = 4;
DD(5) = 400*log(Rstarn);

%% Spreads
ZZ(6,E_Rktil) = 4;
ZZ(6,R_t) = -4;
DD(6) = 400*log(sprd);

%% Output Level
ZZ(7,zlev_t) = 1;
ZZ(7,y_t)    = 1;
DD(7)        = 100*(log(ystar)+gam);

%% Variance matrix of shocks.
QQ(z_sh,z_sh) = sig_z^2;
QQ(phi_sh,phi_sh) = sig_phi^2;
QQ(mu_sh,mu_sh) = sig_mu^2;
QQ(b_sh,b_sh) = sig_b^2;
QQ(g_sh,g_sh) = sig_g^2;
QQ(laf_sh,laf_sh) = sig_laf^2;  
QQ(sigw_sh,sigw_sh) = sig_sigw^2;  
QQ(mue_sh,mue_sh) = sig_mue^2;  
QQ(gamm_sh,gamm_sh) = sig_gamm^2;  
QQ(r_sh,r_sh) = sig_r^2;   %CP: divided by nant - CHANGE

if exist('nant','var')
    if nant > 0
        % These lines set the standard deviations for the anticipated
        % shocks to be equal to the standard deviation for the
        % unanticipated policy shock.
        for i = 1:nant
            eval(strcat('QQ(r_shl',num2str(i),',r_shl',num2str(i),...
                        ') = sig_r_ant(',num2str(i),')^2;'));
        end
    
        for i = 1:nant
            ZZ(7 + i,:) = ZZ(5,:)*(TTT^i);
            DD(7 + i) = 400*log(Rstarn);
        end
    end
end
