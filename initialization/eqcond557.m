%% marginal cost
G0(margcost,mc_t) = 1;
G0(margcost,w_t) = -(1-alp);
G0(margcost,rk_t) = -alp;


%% price setting
G0(prsett,pi_t) = (1+iota_p*bet)*zeta_p;
G0(prsett,E_pi) = -bet*zeta_p;
G0(prsett,mc_t) = -(1-zeta_p)*(1-bet*zeta_p);

  G0(prsett,laf_t) = -1;
  %% transform laf_t, old version: -laf*(1-zeta_p)*(1-bet*zeta_p)/(1+laf);

G1(prsett,pi_t) = iota_p*zeta_p;


%% capital accumulation
G0(capacc,kbar_t) = 1;
G0(capacc,z_t) = (1-istokbarst);  
G0(capacc,i_t) = -istokbarst;
  G0(capacc,mu_t) = -istokbarst*(1+bet)*exp(2*zstar)*s2; 
  %% transformed mu_t; old version -istokbarst

G1(capacc,kbar_t) = (1-istokbarst);


%% effective capital 
G0(effcap,k_t) = 1;
G0(effcap,u_t) = -1;
G0(effcap,z_t) = 1;

G1(effcap,kbar_t) = 1;


%% Euler
G0(euler,xi_t) = (exp(zstar)-h*bet)*(exp(zstar)-h);
    
  G0(euler,b_t) = -(exp(2*zstar)+h^2*bet); %% old version: -(exp(zstar)-h)*exp(zstar);
  G0(euler,b_t) = G0(euler,b_t) + bet*h*rho_b*exp(-zstar)*(exp(2*zstar)+h^2*bet); 
  %% old version: G0(euler,b_t)+(exp(zstar)-h)*h*bet*rho_b; this is Eb_t+1
  
G0(euler,z_t) = h*exp(zstar);
G0(euler,z_t) = G0(euler,z_t)-h*bet*exp(zstar)*rho_z;   %% this is Ez_t+1 
G0(euler,c_t) = (exp(2*zstar)+h^2*bet);
G0(euler,E_c) = -h*bet*exp(zstar);

G1(euler,c_t) = h*exp(zstar);

%% money demand
G0(moneydem,m_t) = nu_m;
G0(moneydem,xi_t) = 1;
G0(moneydem,R_t) = 1/(Rstarn-1);


%% marginal utility
G0(margut,xi_t) = 1;
G0(margut,E_xi) = -1;
G0(margut,z_t) = rho_z;   %% this is Ez_t+1 
G0(margut,R_t) = -1;
G0(margut,E_pi) = 1;


%% capital utilization

  G0(utcap,u_t) = a2;
  G0(utcap,rk_t) = -rkstar;


%% optimal wage
G0(optwage,wtil_t) = 1+nu_l*(1+law)/law;
G0(optwage,w_t) = 1+zeta_w*bet*nu_l*(1+law)/law;
G0(optwage,E_wtil) = -zeta_w*bet*(1+nu_l*(1+law)/law);
G0(optwage,E_w) = -zeta_w*bet*(1+nu_l*(1+law)/law);

  G0(optwage,b_t) = -(1-zeta_w*bet)*(exp(2*zstar)+h^2*bet)*exp(-zstar)/(exp(zstar)-h);
  %% old version: -(1-zeta_w*bet)
  G0(optwage,phi_t) = -1; %% transformed phi_t; old version: -(1-zeta_w*bet);

G0(optwage,L_t) = -(1-zeta_w*bet)*nu_l;
G0(optwage,xi_t) = (1-zeta_w*bet);
G0(optwage,E_pi) = -zeta_w*bet*(1+nu_l*(1+law)/law);

% This line corrected since it seems it should include an iota_w*z_t term.
G0(optwage,z_t) = -zeta_w*bet*rho_z*(1+nu_l*(1+law)/law);
%G0(optwage,z_t) = -zeta_w*bet*(rho_z-iota_w)*(1+nu_l*(1+law)/law);

G0(optwage,pi_t) = zeta_w*bet*iota_w*(1+nu_l*(1+law)/law);


%% aggregate wage evolution
G0(aggwage,w_t) = 1;
G0(aggwage,pi_t) = 1;
G0(aggwage,z_t) = 1;
G0(aggwage,wtil_t) = -(1-zeta_w)/zeta_w;

G1(aggwage,w_t) = 1;
G1(aggwage,pi_t) = iota_w;


%% capital labor ratio
G0(caplabrat,k_t) = 1;
G0(caplabrat,L_t) = -1;
G0(caplabrat,w_t) = -1;
G0(caplabrat,rk_t) = 1;


%% aggregate resources
G0(resources,y_t) = 1;
G0(resources,c_t) = -cstar/(cstar+istar);
G0(resources,i_t) = -istar/(cstar+istar);
G0(resources,g_t) = -1;
G0(resources,u_t) = -rkstar*kstar/(cstar+istar);


%% production function
G0(prod,y_t) = 1;
G0(prod,k_t) = -alp*(ystar+Bigphi)/ystar;
G0(prod,L_t) = -(1-alp)*(ystar+Bigphi)/ystar;


%% Taylor
G0(taylor,R_t) = 1;

G0(taylor,pi_t) = -(1-rho_r)*psi1*1/4;     %% targeting 4-quarter inflation: current quarter
G0(taylor,pi_t1) = -(1-rho_r)*psi1*1/4;    %% targeting 4-quarter inflation: t-1
G0(taylor,pi_t2) = -(1-rho_r)*psi1*1/4;    %% targeting 4-quarter inflation: t-2
G1(taylor,pi_t2) = (1-rho_r)*psi1*1/4;     %% targeting 4-quarter inflation: t-3

%G0(taylor,y_t) = -(1-rho_r)*psi2;     %% targeting output level: current quarter

G0(taylor,y_t) = -(1-rho_r)*psi2*1/4;     %% targeting 4-quarter output growth: current quarter
G1(taylor,y_t) = -(1-rho_r)*psi2*1/4;  
G0(taylor,z_t) = -(1-rho_r)*psi2*1/4;

G1(taylor,gL1_t) = (1-rho_r)*psi2*1/4;  %% targeting 4-quarter output growth: t-1 quarter
G1(taylor,gL2_t) = (1-rho_r)*psi2*1/4;  %% targeting 4-quarter output growth: t-2 quarter
G1(taylor,gL3_t) = (1-rho_r)*psi2*1/4;  %% targeting 4-quarter output growth: t-3 quarter

% Note: Revised "GAM" TO "G"

G1(taylor,R_t) = rho_r;

if exist('nant','var')
    if nant > 0
        G0(taylor,r_t) = -1;
    else
        PSI(taylor,r_sh) = 1;
    end
else
    PSI(taylor,r_sh) = 1;
end

%% GL1
G0(gL1,gL1_t) = 1;
G0(gL1,y_t) = -1;     %% targeting output growth
G1(gL1,y_t) = -1;  
G0(gL1,z_t) = -1;  

%% GL2
G0(gL2,gL2_t) = 1;
G1(gL2,gL1_t) = 1; 

%% GL3
G0(gL3,gL3_t) = 1;
G1(gL3,gL2_t) = 1; 
  
%% pi_t1
G0(pi1,pi_t1) = 1;
G1(pi1,pi_t)  = 1;

%% pi_t2
G0(pi2,pi_t2) = 1;
G1(pi2,pi_t1)  = 1;

%% R_t1
%G0(pi2,R_t1) = 1;
%G1(pi2,R_t)  = 1;


%% investment FOC
G0(invfoc,qk_t) = exp(-2*zstar)/s2;
G0(invfoc,mu_t) = (1+bet); %% old version 1/(s2*exp(2*zstar))
G0(invfoc,z_t) = -1;
G0(invfoc,z_t) = G0(invfoc,z_t)+bet*rho_z;  %% this is Ez_t+1 
G0(invfoc,i_t) = -(1+bet);
G0(invfoc,E_i) = bet;
G1(invfoc,i_t) = -1;


%% return to capital
G0(rettocap,Rktil_t) = 1;
G0(rettocap,pi_t) = -1;
G0(rettocap,rk_t) = -rkstar/(rkstar+1-del);
G0(rettocap,qk_t) = -(1-del)/(rkstar+1-del);
G1(rettocap,qk_t) = -1;


%% spreads
G0(spread,E_Rktil) = 1;
G0(spread,R_t) = -1;
G0(spread,qk_t) = -zeta_spb;
G0(spread,kbar_t) = -zeta_spb;
G0(spread,n_t) = zeta_spb;
G0(spread,sigw_t) = -1;
G0(spread,mue_t) = -1;

%% n evol
G0(nevol,n_t) = 1;
G0(nevol,gamm_t) = -1;
G0(nevol,z_t) = gammstar*vstar/nstar;
G0(nevol,Rktil_t) = -zeta_nRk;
G0(nevol,pi_t) = (zeta_nRk - zeta_nR);
G1(nevol,R_t) = -zeta_nR;
G1(nevol,sigw_t) = -zeta_nsigw/zeta_spsigw;
G1(nevol,mue_t) = -zeta_nmue/zeta_spmue;
G1(nevol,qk_t) = zeta_nqk;
G1(nevol,kbar_t) = zeta_nqk;
G1(nevol,n_t) = zeta_nn;

%% level technology evolution
G0(zlevevol,zlev_t) = 1;
G1(zlevevol,zlev_t) = 1;
G0(zlevevol,z_t)    = -1;
%C(zlevevol) = 100*gam;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% exogenous equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% z shock
G0(n_eqc+z_sh,z_t) = 1;
G1(n_eqc+z_sh,z_t) = rho_z;
PSI(n_eqc+z_sh,z_sh) = 1;

%% phi shock
G0(n_eqc+phi_sh,phi_t) = 1;
G1(n_eqc+phi_sh,phi_t) = rho_phi;
PSI(n_eqc+phi_sh,phi_sh) = 1;


%% mu shock
G0(n_eqc+mu_sh,mu_t) = 1;
G1(n_eqc+mu_sh,mu_t) = rho_mu;
PSI(n_eqc+mu_sh,mu_sh) = 1;


    
 %% b shock
G0(n_eqc+b_sh,b_t) = 1;
G1(n_eqc+b_sh,b_t) = rho_b;
PSI(n_eqc+b_sh,b_sh) = 1;


%% g shock
G0(n_eqc+g_sh,g_t) = 1;
G1(n_eqc+g_sh,g_t) = rho_g;
PSI(n_eqc+g_sh,g_sh) = 1;


%% laf shock
G0(n_eqc+laf_sh,laf_t) = 1;
G1(n_eqc+laf_sh,laf_t) = rho_laf; 
PSI(n_eqc+laf_sh,laf_sh) = 1;


%% sigw shock
G0(n_eqc+sigw_sh,sigw_t) = 1;
G1(n_eqc+sigw_sh,sigw_t) = rho_sigw; 
PSI(n_eqc+sigw_sh,sigw_sh) = 1;


%% mue shock
G0(n_eqc+mue_sh,mue_t) = 1;
G1(n_eqc+mue_sh,mue_t) = rho_mue; 
PSI(n_eqc+mue_sh,mue_sh) = 1;


%% gamm shock
G0(n_eqc+gamm_sh,gamm_t) = 1;
G1(n_eqc+gamm_sh,gamm_t) = rho_gamm; 
PSI(n_eqc+gamm_sh,gamm_sh) = 1;



if exist('nant','var')
    if nant > 0
        %% r shock
        G0(n_eqc+r_sh,r_t) = 1;
        PSI(n_eqc+r_sh,r_sh) = 1;

        % This section adds the anticipated shocks. There is one state for all the
        % anticipated shocks that will hit in a given period (i.e. r_tl2 holds
        % those that will hit in two periods), and the equations are set up so that
        % r_tl2 last period will feed into r_tl1 this period (and so on for other
        % numbers), and last period's r_tl1 will feed into the r_t process (and
        % affect the Taylor Rule this period).

        if nant > 0

            G1(n_eqc+r_sh,r_tl1) = 1;
            G0(n_eqc+r_shl1,r_tl1) = 1;
            PSI(n_eqc+r_shl1,r_shl1) = 1;

            if nant > 1
                for i = 2:nant
                    eval(strcat('G1(n_eqc+r_shl',num2str(i-1),',r_tl',num2str(i),') = 1;'));
                    eval(strcat('G0(n_eqc+r_shl',num2str(i),',r_tl',num2str(i),') = 1;'));
                    eval(strcat('PSI(n_eqc+r_shl',num2str(i),',r_shl',num2str(i),') = 1;'));
                end
            end

        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% expectational equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% E_c
G0(n_eqc+n_exo+Ec_sh,c_t) = 1;
G1(n_eqc+n_exo+Ec_sh,E_c) = 1;
PIE(n_eqc+n_exo+Ec_sh,Ec_sh) = 1;

%% E_pi
G0(n_eqc+n_exo+Epi_sh,pi_t) = 1;
G1(n_eqc+n_exo+Epi_sh,E_pi) = 1;
PIE(n_eqc+n_exo+Epi_sh,Epi_sh) = 1;

%% E_w
G0(n_eqc+n_exo+Ew_sh,w_t) = 1;
G1(n_eqc+n_exo+Ew_sh,E_w) = 1;
PIE(n_eqc+n_exo+Ew_sh,Ew_sh) = 1;

%% E_wtil
G0(n_eqc+n_exo+Ewtil_sh,wtil_t) = 1;
G1(n_eqc+n_exo+Ewtil_sh,E_wtil) = 1;
PIE(n_eqc+n_exo+Ewtil_sh,Ewtil_sh) = 1;

%% E_xi
G0(n_eqc+n_exo+Exi_sh,xi_t) = 1;
G1(n_eqc+n_exo+Exi_sh,E_xi) = 1;
PIE(n_eqc+n_exo+Exi_sh,Exi_sh) = 1;

%% E_i
G0(n_eqc+n_exo+Ei_sh,i_t) = 1;
G1(n_eqc+n_exo+Ei_sh,E_i) = 1;
PIE(n_eqc+n_exo+Ei_sh,Ei_sh) = 1;

%% E_Rktil
G0(n_eqc+n_exo+ERktil_sh,Rktil_t) = 1;
G1(n_eqc+n_exo+ERktil_sh,E_Rktil) = 1;
PIE(n_eqc+n_exo+ERktil_sh,ERktil_sh) = 1;

