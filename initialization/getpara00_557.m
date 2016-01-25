function [alp,zeta_p,iota_p,del,ups,Bigphi,s2,h,a2,nu_l,nu_m,zeta_w,iota_w,law,rstar,psi1,psi2,rho_r,pistar,...
          Fom,sprd,zeta_spb,gammstar,NEW_5,...
          gam,Lstar,chi,laf,gstar,Ladj,...
          rho_z,rho_phi,rho_chi,rho_laf,rho_mu,rho_b,rho_g,rho_sigw,rho_mue,rho_gamm,...
          sig_z,sig_phi,sig_chi,sig_laf,sig_mu,sig_b,sig_g,sig_r,sig_sigw,sig_mue,sig_gamm,...
          bet,zstar,phi,istokbarst,rkstar,cstar,ystar,istar,kstar,kbarstar,Rstarn,wstar,wadj,mstar,...
          zeta_nRk, zeta_nR, zeta_nqk, zeta_nn, zeta_nmue, zeta_spmue, zeta_nsigw, zeta_spsigw, ...
          vstar, nstar,Rkstar, sig_r_ant] = getpara00_557(para) % moved zeta_spd and gammstar higher on the list

nantpad = 20;
      
alp = para(1);
zeta_p = para(2);
iota_p = para(3);
del = para(4);
ups = 1;
Bigphi = para(6);
s2 = para(7);
h = para(8);
a2 = para(9);
nu_l = para(10);
nu_m = para(11);
zeta_w = para(12);
iota_w = para(13);
law = para(14);
bet  = exp(-para(15)/400);  %% (determines bet)
psi1 = para(16);
psi2 = para(17);
rho_r = para(18);
pistar = exp(para(19)/400);
Fom = 1-(1-para(20))^(1/4);  %% F(omega) from annualized to quarterly default prob
sprd = (1+para(21)/100)^(1/4); %exp(para(21)/400);    %% st st spread from annual perc to quarterly number
zeta_spb = para(22);
gammstar = para(23);
NEW_5 = para(24);

npara = 24; 


%% exogenous processes - level

gam = para(npara+1)/400;
wadj = para(npara+2);
Lstar = 1;
chi = para(npara+3);
laf = para(npara+4);
gstar = 1+para(npara+5);
Ladj = para(npara+6);

npara = npara+6;


%% exogenous processes - autocorrelation

rho_z = para(npara+1);
rho_phi = para(npara+2);
rho_chi = para(npara+3);
rho_laf = para(npara+4);
rho_mu = para(npara+5);
rho_b = para(npara+6);
rho_g = para(npara+7);
rho_sigw = para(npara+8);
rho_mue = para(npara+9);
rho_gamm = para(npara+10);

npara = npara+10;    


%% exogenous processes - standard deviation

sig_z = para(npara+1);
sig_phi = para(npara+2);
sig_chi = para(npara+3);
sig_laf = para(npara+4);
sig_mu = para(npara+5);
sig_b = para(npara+6);
sig_g = para(npara+7);
sig_r = para(npara+8);
sig_sigw = para(npara+9);
sig_mue = para(npara+10);
sig_gamm = para(npara+11);

npara = npara+11;    


%% Standard deviations of the anticipated policy shocks
for i = 1:nantpad
eval(strcat('sig_r',num2str(i), '= para(npara +',num2str(i),');'));
eval(strcat('sig_r_ant(', num2str(i),') = sig_r',num2str(i),';'));
end


npara = npara+nantpad;

%% Parameters (implicit) -- from steady state

zstar = gam+(alp\(1-alp))*log(ups); 

rstar = (1/bet)*exp(gam)*ups^(alp/(1-alp)); 

rkstar = sprd*(1/bet)*exp(gam)*ups^(1/(1-alp))-(1-del); %notice that this depends on the spread.

omegastar = (alp^(alp)*(1-alp)^(1-alp)*rkstar^(-alp)/(1+laf))^(1/(1-alp));

%if any(mspec == [105]) 
 %Bigphi = laf*omegastar*Lstar/(1-alp);
 %end

kstar = (alp/(1-alp))*omegastar*Lstar/rkstar;

kbarstar = kstar*exp(gam)*ups^(1/(1-alp));

istokbarst = 1-((1-del)/(exp(gam)*ups^(1/(1-alp))));

istar = kbarstar*istokbarst;

ystar = (kstar^alp)*(Lstar^(1-alp))-Bigphi;
if ystar <= 0

    disp([alp,  bet, kstar,Lstar])
    dm([ystar,Lstar,kstar,Bigphi])
    %keyboard;
end
cstar = (ystar/gstar)-istar;

xistar = (1/cstar)*((1/(1-h*exp(-gam)*ups^(-alp/(1-alp))))-(h*bet/(exp(gam)*ups^(alp/(1-alp))-h)));

phi = Lstar^(-nu_l)*omegastar*xistar/(1+law);

Rstarn = pistar*rstar;

wstar = ( 1/(1+laf) * (alp^alp) * (1-alp)^(1-alp) * rkstar^(-alp) )^(1/(1-alp));

mstar = ( chi * Rstarn/(Rstarn-1) * xistar^(-1) )^(1 / nu_m);

%% Solve for financial frictions parameters

% Fom = 0.03; Fom = 1-(1-Fom)^(1/4);
% sprd = 2; sprd = (1+sprd/100)^(1/4);
% zeta_spb = 0.045;
% gammstar = 0.99;

% zwstar = norminv(Fom);
% sigwstar = (0.28/4)^(1/2);
% muestar = mufcn(zwstar,sigwstar,sprd)
% zetaspb = zetaspbfcn(zwstar,sigwstar,sprd)
% nkstar = nkfcn(zwstar,sigwstar,sprd)
% knstar = 1/nkstar
% % sig = 0.1:.1:1;
% % for j=1:length(sig)
% %     zetaspb(j) = zetaspbfcn(zwstar,sig(j),sprd);
% % end
% % plot(sig,zetaspb)
% return


% solve for sigmaomegastar and zomegastar
zwstar = norminv(Fom);
sigwstar = fzero(@(sigma)zetaspbfcn(zwstar,sigma,sprd)-zeta_spb,0.5);
% zetaspbfcn(zwstar,sigwstar,sprd)-zeta_spb % check solution

% evaluate omegabarstar
omegabarstar = omegafcn(zwstar,sigwstar);

% evaluate all BGG function elasticities
Gstar = Gfcn(zwstar,sigwstar);
Gammastar = Gammafcn(zwstar,sigwstar);
dGdomegastar = dGdomegafcn(zwstar,sigwstar);
d2Gdomega2star = d2Gdomega2fcn(zwstar,sigwstar);
dGammadomegastar = dGammadomegafcn(zwstar);
d2Gammadomega2star = d2Gammadomega2fcn(zwstar,sigwstar);
dGdsigmastar = dGdsigmafcn(zwstar,sigwstar);
d2Gdomegadsigmastar = d2Gdomegadsigmafcn(zwstar,sigwstar);
dGammadsigmastar = dGammadsigmafcn(zwstar,sigwstar);
d2Gammadomegadsigmastar = d2Gammadomegadsigmafcn(zwstar,sigwstar);

% evaluate mu, nk, and Rhostar
muestar = mufcn(zwstar,sigwstar,sprd);
nkstar = nkfcn(zwstar,sigwstar,sprd);
Rhostar = 1/nkstar-1;

% evaluate wekstar and vkstar
wekstar = (1-gammstar/bet)*nkstar...
    -gammstar/bet*(sprd*(1-muestar*Gstar)-1);
vkstar = (nkstar-wekstar)/gammstar;

% evaluate nstar and vstar
nstar = nkstar*kstar;
vstar = vkstar*kstar;

% a couple of combinations
GammamuG = Gammastar-muestar*Gstar;
GammamuGprime = dGammadomegastar-muestar*dGdomegastar;

% elasticities wrt omegabar
zeta_bw = zetabomegafcn(zwstar,sigwstar,sprd);
zeta_zw = zetazomegafcn(zwstar,sigwstar,sprd);
zeta_bw_zw = zeta_bw/zeta_zw;

% elasticities wrt sigw
zeta_bsigw = sigwstar*(((1-muestar*dGdsigmastar/dGammadsigmastar)/...
    (1-muestar*dGdomegastar/dGammadomegastar)-1)*dGammadsigmastar*sprd+...
    muestar*nkstar*(dGdomegastar*d2Gammadomegadsigmastar-dGammadomegastar*d2Gdomegadsigmastar)/...
    GammamuGprime^2)/...
    ((1-Gammastar)*sprd+dGammadomegastar/GammamuGprime*(1-nkstar));
zeta_zsigw = sigwstar*(dGammadsigmastar-muestar*dGdsigmastar)/GammamuG;
zeta_spsigw = (zeta_bw_zw*zeta_zsigw-zeta_bsigw)/(1-zeta_bw_zw);

% elasticities wrt mue
zeta_bmue = muestar*(nkstar*dGammadomegastar*dGdomegastar/GammamuGprime+dGammadomegastar*Gstar*sprd)/...
    ((1-Gammastar)*GammamuGprime*sprd+dGammadomegastar*(1-nkstar));
zeta_zmue = -muestar*Gstar/GammamuG;
zeta_spmue = (zeta_bw_zw*zeta_zmue-zeta_bmue)/(1-zeta_bw_zw);

% some ratios/elasticities
Rkstar = sprd*pistar*rstar; % (rkstar+1-delta)/ups*pistar;
zeta_Gw = dGdomegastar/Gstar*omegabarstar;
zeta_Gsigw = dGdsigmastar/Gstar*sigwstar;

% % elasticities for equity equation
% zeta_vRk = Rkstar/pistar/exp(zstar)/vkstar*(1-muestar*Gstar*(1-zeta_Gw/zeta_zw));
% zeta_vR = 1/bet/vkstar*(1-nkstar+muestar*Gstar*sprd*zeta_Gw/zeta_zw);
% zeta_vqk = Rkstar/pistar/exp(zstar)/vkstar*(1-muestar*Gstar*...
%     (1+zeta_Gw/zeta_zw*nkstar/(1-nkstar)))-1/bet/vkstar;
% zeta_vn = 1/bet*nkstar/vkstar+...
%     Rkstar/pistar/exp(zstar)/vkstar*muestar*Gstar*zeta_Gw/zeta_zw/Rhostar;
% zeta_vmue = Rkstar/pistar/exp(zstar)/vkstar*muestar*Gstar*(1-zeta_Gw*zeta_zmue/zeta_zw);
% zeta_vsigw = Rkstar/pistar/exp(zstar)/vkstar*muestar*Gstar*(zeta_Gsigw-zeta_Gw/zeta_zw*zeta_zsigw);

% elasticities for the net worth evolution
zeta_nRk = gammstar*Rkstar/pistar/exp(zstar)*(1+Rhostar)*(1-muestar*Gstar*(1-zeta_Gw/zeta_zw));
zeta_nR = gammstar/bet*(1+Rhostar)*(1-nkstar+muestar*Gstar*sprd*zeta_Gw/zeta_zw);
zeta_nqk = gammstar*Rkstar/pistar/exp(zstar)*(1+Rhostar)*(1-muestar*Gstar*(1+zeta_Gw/zeta_zw/Rhostar))...
    -gammstar/bet*(1+Rhostar);
zeta_nn = gammstar/bet+gammstar*Rkstar/pistar/exp(zstar)*(1+Rhostar)*muestar*Gstar*zeta_Gw/zeta_zw/Rhostar;
zeta_nmue = gammstar*Rkstar/pistar/exp(zstar)*(1+Rhostar)*muestar*Gstar*(1-zeta_Gw*zeta_zmue/zeta_zw);
zeta_nsigw = gammstar*Rkstar/pistar/exp(zstar)*(1+Rhostar)*muestar*Gstar*(zeta_Gsigw-zeta_Gw/zeta_zw*zeta_zsigw);

% % show some ratios...
% v.Fom = Fom;
% v.sprd = sprd;
% v.zeta_spb = zeta_spb;
% v.gammstar = gammstar;
% 
% v.sigwstar = sigwstar;
% v.muestar = muestar;
% v.nkstar = nkstar;
% v.vkstar = vkstar;
% v.wekstar = wekstar;
% v.Rhostar = Rhostar;
% 
% v.zeta_nRk = zeta_nRk;
% v.zeta_nR = zeta_nR;
% v.zeta_nqk = zeta_nqk;
% v.zeta_nn = zeta_nn;
% v.zeta_nmue = zeta_nmue;
% v.zeta_spmue = zeta_spmue;
% v.zeta_nsigw = zeta_nsigw;
% v.zeta_spsigw = zeta_spsigw;
% v.vstar = vstar;
% v.nstar = nstar;
% 
% v.zwstar = zwstar;
% v.omegabarstar = omegabarstar;
% v.Gstar = Gstar;
% v.Gammastar = Gammastar;
% v.dGdomegastar = dGdomegastar;
% v.d2Gdomega2star = d2Gdomega2star;
% v.dGammadomegastar = dGammadomegastar;
% v.d2Gammadomega2star = d2Gammadomega2star;
% v.dGdsigmastar = dGdsigmastar;
% v.d2Gdomegadsigmastar = d2Gdomegadsigmastar;
% v.dGammadsigmastar = dGammadsigmastar;
% v.d2Gammadomegadsigmastar = d2Gammadomegadsigmastar;
% 
% v


end

function f=zetaspbfcn(z,sigma,sprd)
zetaratio = zetabomegafcn(z,sigma,sprd)/zetazomegafcn(z,sigma,sprd);
nk = nkfcn(z,sigma,sprd);
f = -zetaratio/(1-zetaratio)*nk/(1-nk);
end

function f=zetabomegafcn(z,sigma,sprd)
nk = nkfcn(z,sigma,sprd);
mustar = mufcn(z,sigma,sprd);
omegastar = omegafcn(z,sigma);
Gammastar = Gammafcn(z,sigma);
Gstar = Gfcn(z,sigma);
dGammadomegastar = dGammadomegafcn(z);
dGdomegastar = dGdomegafcn(z,sigma);
d2Gammadomega2star = d2Gammadomega2fcn(z,sigma);
d2Gdomega2star = d2Gdomega2fcn(z,sigma);
f = omegastar*mustar*nk*(d2Gammadomega2star*dGdomegastar-d2Gdomega2star*dGammadomegastar)/...
    (dGammadomegastar-mustar*dGdomegastar)^2/sprd/...
    (1-Gammastar+dGammadomegastar*(Gammastar-mustar*Gstar)/(dGammadomegastar-mustar*dGdomegastar));
end

function f=zetazomegafcn(z,sigma,sprd)
mustar = mufcn(z,sigma,sprd);
f = omegafcn(z,sigma)*(dGammadomegafcn(z)-mustar*dGdomegafcn(z,sigma))/...
    (Gammafcn(z,sigma)-mustar*Gfcn(z,sigma));
end

function f=nkfcn(z,sigma,sprd)
f = 1-(Gammafcn(z,sigma)-mufcn(z,sigma,sprd)*Gfcn(z,sigma))*sprd;
end

function f=mufcn(z,sigma,sprd)
f = (1-1/sprd)/(dGdomegafcn(z,sigma)/dGammadomegafcn(z)*(1-Gammafcn(z,sigma))+Gfcn(z,sigma));
end

function f=omegafcn(z,sigma)
f = exp(sigma*z-1/2*sigma^2);
end

function f=Gfcn(z,sigma)
f = normcdf(z-sigma);
end

function f=Gammafcn(z,sigma)
f = omegafcn(z,sigma)*(1-normcdf(z))+normcdf(z-sigma);
end

function f=dGdomegafcn(z,sigma)
f=normpdf(z)/sigma;
end

function f=d2Gdomega2fcn(z,sigma)
f = -z*normpdf(z)/omegafcn(z,sigma)/sigma^2;
end

function f=dGammadomegafcn(z)
f = 1-normcdf(z);
end

function f=d2Gammadomega2fcn(z,sigma)
f = -normpdf(z)/omegafcn(z,sigma)/sigma;
end

function f=dGdsigmafcn(z,sigma)
f = -z*normpdf(z-sigma)/sigma;
end

function f=d2Gdomegadsigmafcn(z,sigma)
f = -normpdf(z)*(1-z*(z-sigma))/sigma^2;
end

function f=dGammadsigmafcn(z,sigma)
f = -normpdf(z-sigma);
end

function f=d2Gammadomegadsigmafcn(z,sigma)
f = (z/sigma-1)*normpdf(z);
end

