function [para,para_names,para_mask,para_fix,npara,polipar,polivalue,bounds] = mspec_parameters_557(subspec,dataset)

%% names of parameters
% the parameters in para_names must be lsited in the order they appear below and in the model's priorfile

para_names = { '\alpha';'\zeta_p';'\iota_p';'\delta';'\Upsilon';'\Phi';'S''''';'h';'a''''';'\nu_l';'\nu_m';'\zeta_w';...
               '\iota_w';'\lambda_w';'r^*';...
               '\psi_1';'\psi_2';'\rho_{r}';'\pi^*';'F(\omega)';'spr_*';'\zeta_{sp}';'\gamma_*';'NEW_5';...
               '\gamma';'W^{adj}';'\chi';'\lambda_f';'g^*';'L^{adj}';...
               '\rho_{z}';'\rho_{\phi}';'\rho_{\chi}';'\rho_{\lambda_f}';'\rho_{\mu}';'\rho_{b}';'\rho_{g}';'\rho_{sigw}';'\rho_{mue}';'\rho_{gamm}';...
               '\sigma_{z}';'\sigma_{\phi}';'\sigma_{\chi}';'\sigma_{\lambda_f}';'\sigma_{\mu}';'\sigma_{b}';'\sigma_{g}';'\sigma_{r}';...
               '\sigma_{sigw}';'\sigma_{mue}';'\sigma_{gamm}'};

%% starting values for parameters
para = zeros(100,1);
nantpad = 20;

% para_mask is npara x 1
%   to fix a parameter under specification mspec, set the corresponding element of para_mask to 1
% para_fix is npara x 1
%   and contains the values for the fixed parameters (all other entries are irrelevant)

para_mask = zeros(100,1);
para_fix  = zeros(100,1);

para(1) = .33;     %para_mask(1) = 1;  para_fix(1) = para(1);              %% alp;         1
para(2) = .75;     %para_mask(2) = 1;  para_fix(2) = para(2);              %% zeta_p;      2
para(3) = .5;      %para_mask(3) = 1;   para_fix(3) = para(3);              %% iota_p;      3  
para(4) = .025;     para_mask(4) = 1;   para_fix(4) = para(4);              %% del;         4
para(5) = 963;        %para_mask(5) = 1;   para_fix(5) = para(5);              %% Ztil_0; .1      5
para(6) = 0;        para_mask(6) = 1;   para_fix(6) = para(6);              %% Bigphi;      6
para(7) = 4;       %para_mask(7) = 1;  para_fix(7) = para(7);              %% s2;           7
para(8) = .7;      %para_mask(8) = 1;  para_fix(8) = para(8);              %% h;           8
para(9) = .2;      %para_mask(9) = 1;  para_fix(9) = para(9);              %% a2;          9
para(10) = 2;      %para_mask(10) = 1;  para_fix(10) = para(10);           %% nu_l         10
para(11) = 2;       para_mask(11) = 1;  para_fix(11) = para(11);            %% nu_m         11
para(12) = .8;     %para_mask(12) = 1; para_fix(12) = para(12);            %% zeta_w;      12
para(13) = .5;    %para_mask(13) = 1;   para_fix(13) = para(13);          %% iota_w;      13
para(14) = .3;      para_mask(14) = 1;  para_fix(14) = para(14);            %% law;         14
para(15) = 2;      %para_mask(15) = 1; para_fix(15) = para(15);            %% beta;        15
para(16) = 1.7;    %para_mask(16) = 1; para_fix(16) = para(16);            %% psi1;        16
para(17) = .125;   %para_mask(17) = 1; para_fix(17) = para(17);            %% psi2;        17
para(18) = .8;     %para_mask(18) = 1; para_fix(18) = para(18);            %% rho_r;       18
para(19) = 3;      %para_mask(19) = 1; para_fix(19) = para(19);            %% pistar;      19

npara = 19;

para(npara+1) = .03; para_mask(20) = 1;  para_fix(20) = para(20);              %% Fom
para(npara+2) = 2; %para_mask(21) = 1;  para_fix(21) = para(21);              %% st st spread
para(npara+3) = 0.05; %para_mask(22) = 1;  para_fix(22) = para(22);              %% zeta_sp
para(npara+4) = .99; para_mask(23) = 1;  para_fix(23) = para(23);              %% gammstar
para(npara+5) = 0; para_mask(24) = 1;  para_fix(24) = para(24);              %% NEW_5
% para(npara+2) = 2; para_mask(21) = 1;  para_fix(21) = para(21);              %% st st spread
% para(npara+3) = .05; %para_mask(22) = 1;  para_fix(22) = para(22);              %% zeta_sp

npara = npara+5;
    
%% exogenous processes - level
if subspec == 551
  para(npara+1) = 1.5;      para_mask(npara+1) = 1;    para_fix(npara+1) = para(npara+1);      %% gam;
elseif subspec == 552
  para(npara+1) = 2;      para_mask(npara+1) = 1;    para_fix(npara+1) = para(npara+1);      %% gam;
else
  para(npara+1) = 2;      %para_mask(npara+1) = 1;    para_fix(npara+1) = para(npara+1);      %% gam;
end
para(npara+2) =  0;      para_mask(npara+2) = 1;    para_fix(npara+2) = para(npara+2);      %% Wadj;
para(npara+3) = .1;      para_mask(npara+3) = 1;    para_fix(npara+3) = para(npara+3);      %% chi;     
para(npara+4) = .15;     para_mask(npara+4) = 1;    para_fix(npara+4) = para(npara+4);      %% laf;     
para(npara+5) = .15;    %para_mask(npara+5) = 1;    para_fix(npara+5) = para(npara+5);      %% gstar;
para(npara+6) = 252;    %para_mask(npara+6) = 1;    para_fix(npara+6) = para(npara+6);      %% Ladj

npara = npara+6;    

%% exogenous processes - autocorrelation
  
para(npara+1) = .95;  %para_mask(npara+1) = 1;    para_fix(npara+1) = para(npara+1);      %% rho_z (rho_zP if mspec == 8);
para(npara+2) = .4; %para_mask(npara+2) = 1;    para_fix(npara+2) = para(npara+2);      %% rho_phi;
para(npara+3) = .9; para_mask(npara+3) = 1;     para_fix(npara+3) = para(npara+3);      %% rho_chi;
para(npara+4) = .9; %para_mask(npara+4) = 1;    para_fix(npara+4) = para(npara+4);      %% rho_laf;
para(npara+5) = .5; %para_mask(npara+5) = 1;    para_fix(npara+5) = para(npara+5);      %% rho_mu;
para(npara+6) = .8; para_mask(npara+6) = 1;     para_fix(npara+6) = para(npara+6);      %% rho_b;
para(npara+7) = .9; %para_mask(npara+7) = 1;    para_fix(npara+7) = para(npara+7);      %% rho_g;

npara = npara+7;

para(npara+1) = .75;  %para_mask(npara+1) = 1;     para_fix(npara+1) = para(npara+1);      %% rho_sigw
para(npara+2) = .75;  para_mask(npara+2) = 1;     para_fix(npara+2) = para(npara+2);      %% rho_mue
para(npara+3) = .75;  para_mask(npara+3) = 1;     para_fix(npara+3) = para(npara+3);      %% rho_gamm

npara = npara+3;

para(npara+1) = .4;     %para_mask(npara+1) = 1;    para_fix(npara+1) = para(npara+1);      %% sig_z (sig_zP if mspec == 8);
para(npara+2) =  1; %para_mask(npara+2) = 1;    para_fix(npara+2) = para(npara+2);      %% sig_phi;
para(npara+3) = 0;  para_mask(npara+3) = 1;     para_fix(npara+3) = para(npara+3);      %% sig_chi;
para(npara+4) = 1;  %para_mask(npara+4) = 1;    para_fix(npara+4) = para(npara+4);      %% sig_laf;
para(npara+5) = 1;  %para_mask(npara+5) = 1;    para_fix(npara+5) = para(npara+5);      %% sig_mu;
para(npara+6) = 0;    para_mask(npara+6) = 1;    para_fix(npara+6) = para(npara+6);  %% sig_b;
para(npara+7) = .30;    %para_mask(npara+7) = 1;        para_fix(npara+7) = para(npara+7);  %% sig_g;
para(npara+8) = .10;    %para_mask(npara+8) = 1;        para_fix(npara+8) = para(npara+8);  %% sig_r;

para(npara+9) = .20/4;  %para_mask(npara+9) = 1;     para_fix(npara+9) = para(npara+9);      %% sig_sigw
para(npara+10) = 0;  para_mask(npara+10) = 1;     para_fix(npara+10) = para(npara+10);      %% sig_mue
para(npara+11) = 0;  para_mask(npara+11) = 1;     para_fix(npara+11) = para(npara+11);      %% sig_gamm

npara = npara+11;

% Standard Deviations of the anticipated policy shocks

for i = 1:nantpad
    eval(strcat('para(npara +',num2str(i),') = 0.20;'));
    if i >=13
    eval(strcat('para_mask(npara +',num2str(i),') = 1;'));
    eval(strcat('para_fix(npara +',num2str(i),') = 0;'));
    end
end

npara = npara+nantpad;

para(3) = 0;        para_mask(3) = 1;    para_fix(3) = para(3);              %% iota_p;      3
para(13) = 0;       para_mask(13) = 1;   para_fix(13) = para(13);            %% iota_w;      13

para = para(1:npara);
para_mask = para_mask(1:npara);
para_fix = para_fix(1:npara);

%% number of parameters in DSGE model

para = para.*(1-para_mask)+para_fix.*para_mask;

%% identify policy and non policy parameters
polipar = 16;
polivalue = 2;

%% bounds for MH

bounds = zeros(100,2);

bounds(1,:) = [1E-5 .99999];    %% alp;         1
bounds(2,:) = [1E-5 .99999];    %% zeta_p;      2
bounds(3,:) = [0 1];            %% iota_p;      3
bounds(4,:) = [1E-5 .99999];    %% del;         4
bounds(5,:) = [0    1E+5];     %% ups;            5
bounds(6,:) = [0    5];     %% Bigphi;          6
bounds(7,:) = [0    20];        %% s2;          7
bounds(8,:) = [1E-5 .99999];    %% h;           8
bounds(9,:) = [0    .99999];    %% a2;          9
bounds(10,:) = [1E-5 10];        %% nu_l        10
bounds(11,:) = [1E-5    100];        %% nu_m    11
bounds(12,:) = [1E-5    .99999];    %% zeta_w;  12
bounds(13,:) = [0 1];            %% iota_w;     13
bounds(14,:) = [1E-5    50];        %% law;     14
bounds(15,:) = [1E-5    10];        %% beta;    15
bounds(16,:) = [1E-5    10];        %% psi1;    16
bounds(17,:) = [1E-5    10];        %% psi2;    17
bounds(18,:) = [1E-5    .99999];    %% rho_r;   18
bounds(19,:) = [-1  10];        %% pistar;      19

npara = 19;

bounds(npara+1,:) = [1E-5   .99999]; %% Fom
bounds(npara+2,:) = [0 100]; %% st st spread
bounds(npara+3,:) = [1E-5   .99999]; %% zeta_sp
bounds(npara+4,:) = [1E-5   .99999]; %% gammstar
bounds(npara+5,:) = [-10000 10000]; %% NEW_5

npara = npara+5;

%% exogenous processes - level
bounds(npara+1,:) = [1E-6   10];        %% gam;     18
bounds(npara+2,:) = [0      10];        %% Lstar;   19
bounds(npara+3,:) = [1E-6   10];        %% chi;
bounds(npara+4,:) = [1E-5   50];        %% laf;     21
bounds(npara+5,:) = [1E-5   .99999];    %% gstar;   22
bounds(npara+6,:) = [1E-6   5000];        %% Ladj     23

npara = npara+6;

%% exogenous processes - autocorrelation
bounds(npara+1,:) = [0      .99999];    %% rho_z (rho_zP if mspec == 8);   24
bounds(npara+2,:) = [1E-5   .99999];    %% rho_phi;
bounds(npara+3,:) = [1E-5   .99999];    %% rho_chi;
bounds(npara+4,:) = [1E-5   .99999];    %% rho_laf;
bounds(npara+5,:) = [1E-5   .99999];    %% rho_mu;
bounds(npara+6,:) = [1E-5   .99999];    %% rho_b;
bounds(npara+7,:) = [1E-5   .99999];    %% rho_g;

npara = npara+7;

bounds(npara+1,:) = [1E-5   .99999]; %% rho_sigw
bounds(npara+2,:) = [1E-5   .99999]; %% rho_mue
bounds(npara+3,:) = [1E-5   .99999]; %% rho_gamm

npara = npara+3;

%% exogenous processes - standard deviation
bounds(npara+1,:) = [1E-7   100];       %% sig_z (sig_zP if mspec == 8);
bounds(npara+2,:) = [1E-7   100];       %% sig_phi;
bounds(npara+3,:) = [1E-7   100];       %% sig_chi;
bounds(npara+4,:) = [1E-7   10000];     %% sig_laf;
bounds(npara+5,:) = [1E-7   100];       %% sig_mu;
bounds(npara+6,:) = [1E-7   100];       %% sig_b;
bounds(npara+7,:) = [1E-7   100];       %% sig_g;
bounds(npara+8,:) = [1E-7   100];       %% sig_r;

npara = npara+8;

bounds(npara+1,:) = [1E-7   100]; %% sig_sigw
bounds(npara+2,:) = [1E-7   100]; %% sig_mue
bounds(npara+3,:) = [1E-7   100]; %% sig_gamm

npara = npara+3;

%% Standard Deviations of the Anticipated Shocks
for i = 1:nantpad
    eval(strcat('bounds(npara +',num2str(i),',:) = [1E-7   100];'));
end

npara = npara+nantpad;

bounds = bounds(1:npara,:);
