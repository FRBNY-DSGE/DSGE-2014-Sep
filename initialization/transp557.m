%% Transformations:
%% format: type, a, b, c
%% Type 1: (transformation for normal)
%% x is [a,b] -> [-1,1] -> [-inf,inf] by (1/c)*c*z/sqrt(1-c*z^2)
%% Type 2:
%% x is [0,inf] -> [-inf,inf] by b + (1/c)*ln(para[i]-a);

function trspec = transp557()

    nantpad = 20;
	trspec = zeros(100,4);
	
	
	trspec(1,:) = [1	1E-5	.99	1]; 	%% alp;
	trspec(2,:) = [1	1E-5	.99	1];   	%% zeta_p;
	trspec(3,:) = [1	1E-5	.99	1];	%% iota_p;
	trspec(4,:) = [1	1E-5	.99	1]; 	%% del;
	trspec(5,:) = [0    0       0   0];	%% Ztil_0;
	trspec(6,:) = [2	1E-5	0	1];	%% Bigphi;
	trspec(7,:) = [2	1E-5	0	1];	%% s2;
	trspec(8,:) = [1	1E-5	.99	1]; 	%% h;
	trspec(9,:) = [2	1E-5	0	1];	%% a2;
	trspec(10,:) = [2	1E-5	0	1];	%% nu_l;
	trspec(11,:) = [2	1E-5	0	1];	%% nu_m;
	trspec(12,:) = [1	1E-5	.99	1];	%% zeta_w;
	trspec(13,:) = [1	1E-5	.99	1];	%% iota_p;
	trspec(14,:) = [2	1E-5	0	1];	%% law;
	trspec(15,:) = [2	1E-5	0	1];	%% rstar;
	trspec(16,:) = [2	1E-5	0	1];	%% psi1;
	trspec(17,:) = [2	1E-5	0	1];	%% psi2;
	trspec(18,:) = [1	1E-5	.99	1];	%% rho_r;
	trspec(19,:) = [0	0	0	0];	%% pistar;
	trspec(20,:) = [1	1E-5	.99	1];     %% F(omega)
	trspec(21,:) = [2	1E-5	0	1];     %% st st spread
    trspec(22,:) = [1	1E-5	.99	1];     %% zeta_sp
    trspec(23,:) = [1	1E-5	.99	1];     %% gammst
    trspec(24,:) = [0       0       0       0]; %% NEW_5
    
	npara = 24;
	
	%% exogenous processes - level
	trspec(npara+1,:) = [2	1E-5	0	1];	%% gam;
	trspec(npara+2,:) = [0	0	0	0];	%% Lstar;
	trspec(npara+3,:) = [2	1E-5	0	1];	%% chi;
	trspec(npara+4,:) = [2	1E-5	0	1];	%% laf;
	trspec(npara+5,:) = [1	1E-5	.99	1];	%% gstar;
	trspec(npara+6,:) = [0	0	0	0];	%% Ladj
	
	npara = npara+6;	
	
	%% exogenous processes - autocorrelation
	trspec(npara+1,:) = [1	1E-5	.99	1];	%% rho_z;
	trspec(npara+2,:) = [1	1E-5	.99	1];	%% rho_phi;
	trspec(npara+3,:) = [1	1E-5	.99	1];	%% rho_chi;
	trspec(npara+4,:) = [1	1E-5	.99	1];	%% rho_laf;
	trspec(npara+5,:) = [1	1E-5	.99	1];	%% rho_mu;
	trspec(npara+6,:) = [1	1E-5	.99	1];	%% rho_b;
	trspec(npara+7,:) = [1	1E-5	.99	1];	%% rho_g;
    trspec(npara+8,:) = [1	1E-5	.99	1];	%% rho_sigw;
    trspec(npara+9,:) = [1	1E-5	.99	1];	%% rho_mue;
    trspec(npara+10,:) = [1	1E-5	.99	1];	%% rho_gamm;
	
	npara = npara+10;	
	
	%% exogenous processes - standard deviation
	trspec(npara+1,:) = [2	1E-5	0	1];	%% sig_z;
	trspec(npara+2,:) = [2	1E-5	0	1];	%% sig_phi;
	trspec(npara+3,:) = [2	1E-5	0	1];	%% sig_chi;
	trspec(npara+4,:) = [2	1E-5	0	1];	%% sig_laf;
	trspec(npara+5,:) = [2	1E-5	0	1];	%% sig_mu;
	trspec(npara+6,:) = [2	1E-5	0	1];	%% sig_b;
	trspec(npara+7,:) = [2	1E-5	0	1];	%% sig_g;
	trspec(npara+8,:) = [2	1E-5	0	1];	%% sig_r;
    trspec(npara+9,:) = [2	1E-5	0	1];	%% sig_sigw;
    trspec(npara+10,:) = [2	1E-5	0	1];	%% sig_mue;
    trspec(npara+11,:) = [2	1E-5	0	1];	%% sig_gamm;
    
    npara = npara+11;
    
    %% Standard Deviation of the Anticipated Shocks
    
    for i = 1:nantpad
        eval(strcat('trspec(npara +',num2str(i),',:) = [2 1E-5  0  1];'));
    end
    
    npara = npara+nantpad;
    
	trspec = trspec(1:npara,:);



