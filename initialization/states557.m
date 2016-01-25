%% Endogenous variables 
y_t = 1;
c_t = 2;
i_t = 3;
k_t = 4;
kbar_t = 5;
L_t = 6;
m_t = 7;
mc_t = 8;
pi_t = 9;
R_t = 10;
rk_t = 11;
u_t = 12;
w_t = 13;
wtil_t = 14;
xi_t = 15;
qk_t = 16;
gL1_t = 17;
gL2_t = 18;
gL3_t = 19;
pi_t1 = 20;
pi_t2 = 21;
Rktil_t = 22;
n_t = 23;
zlev_t  = 24;
n_end = 24;


%% Exogenous variables (equations) and exogenous shocks (if shocks are
%% iid the # of exog equations is less than the # of shocks)
%% note what we call z_t here is really z^*_t in the paper
z_sh = 1;
z_t = n_end+z_sh;
phi_sh = 2;
phi_t = n_end+phi_sh;
mu_sh = 3;
mu_t = n_end+mu_sh;
b_sh = 4;
b_t = n_end+b_sh;
g_sh = 5;
g_t = n_end+g_sh;
laf_sh = 6; %% this is actually iid: for the time being I treat it as an eq.- in case we want to change it
laf_t = n_end+laf_sh;
sigw_sh = 7; %% this is actually iid: for the time being I treat it as an eq.- in case we want to change it
sigw_t = n_end+sigw_sh;
mue_sh = 8; %% this is actually iid: for the time being I treat it as an eq.- in case we want to change it
mue_t = n_end+mue_sh;
gamm_sh = 9; %% this is actually iid: for the time being I treat it as an eq.- in case we want to change it
gamm_t = n_end+gamm_sh;

r_sh = 10;   %% this is the taylor rule shock - iid

if exist('nant','var')
  if nant > 0

    % Set up an r_t state so it can catch last period's one-period-ahead anticipated
    % shocks.
    r_t = n_end + r_sh;

    % These are the anticipated shocks. For each there is both an innovation
    % (for new anticipated shocks, calculated in period T only),
    % and a process, so that the shocks can be passed from period to
    % period.

    for i = 1:nant
      eval(strcat('r_shl',num2str(i),' = ',num2str(r_sh + i),';'));
      eval(strcat('r_tl',num2str(i),' = n_end + r_shl',num2str(i),';'));
    end

    n_exo = 10 + nant;
    nex = 10 + nant; 

  else

    n_exo = 9;
    nex = 10; 

  end  
  
else
  
  n_exo = 9;
  nex = 10;
  
end

%% Expectation terms 
Ec_sh = 1; 
E_c = n_end+n_exo+Ec_sh;
Epi_sh = 2;
E_pi = n_end+n_exo+Epi_sh;
Ew_sh = 3;
E_w = n_end+n_exo+Ew_sh;
Ewtil_sh = 4;
E_wtil = n_end+n_exo+Ewtil_sh;
Exi_sh = 5;
E_xi = n_end+n_exo+Exi_sh;
Ei_sh = 6;
E_i = n_end+n_exo+Ei_sh;
ERktil_sh = 7;
E_Rktil = n_end+n_exo+ERktil_sh;

n_exp = 7;
%% number of endogenous shocks
nend = n_exp;

%% Number of States
nstates = nend + n_end;
