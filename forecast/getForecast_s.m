function [yypred,yypred_s] = getForecast_s(qahead,nvar,nshocks,CCC,TTT,RRR,DD,ZZ,QQ,z,...
                                  ind_r,ind_r_sh,zerobound,bdd_int_rate,nant,sflag,A,nplotstates,varargin)
                              
%% Check inputs      
switch nargin
    case 18
        tflag = 0;
        df = [];
    case 19
        tflag = 0;
        df = [];
        mspec = varargin{1};
    case 20
        tflag = varargin{1}; 
        df = varargin{2};    
    case 21
        tflag = varargin{1};
        df = varargin{2};
        mspec = varargin{3};
    otherwise 
        fprintf('Incorrect number of inputs for getForecast.m');
end

if bdd_int_rate && ( isempty(ind_r) || isempty(ind_r_sh) )
    error('accounting for 0.25 interest rate bound: interest rate not in observables or monetary policy shock not in shocks')
end


%% Initialize output

yypred = zeros(qahead,nvar);
yypred_s = zeros(qahead,nplotstates);

%% Draw shocks for forecasting
if tflag 
    htil_draw = chi2rnd(repmat(df,1,qahead))./repmat(df,1,qahead);
    Shocks = randn(qahead,nshocks)*sqrt(QQ)./sqrt(htil_draw');
else
    Shocks = repmat(sqrt(diag(QQ)'),qahead,1).*randn(qahead,nshocks);
    Shocks(:,ind_r_sh+1:ind_r_sh+nant) = zeros(qahead,nant);
end

%TURN SHOCKS OFF
if sflag
    Shocks = zeros(size(Shocks));
end



%% Calculate forecast 
% Applying state transition equation and measurement equation 
% S_t = TTT*S_(t-1) + RRR*eps_t
% y_t = ZZ*S_t + DD


for t = 1:qahead;
    z_test = CCC+TTT*z+RRR*Shocks(t,:)';
    yypred(t,:) = (DD+ZZ*z_test)';
    yypred_s(t,:) = (A*z_test)';
%    bdd_int_rate = 1; % CHANGE THIS
    

    % This changes the monetary policy shock to account for the 0.25 interest rate bound    
    % For zero bound forecast we need forecasts with bounded interest rate on and off
    
    %THIS IS A HACKKKK (RH: 2012-07-02)
    if ~exist('mspec','var') || ~any(mspec==[803 805 904 828 825 826 827 829 830 835 8351 8352 836 8361 8362 906 8281:8285 8286 8287])
        %ZeroBound=0.25;
        ZeroBound=0.11;
        %ZeroBoundTest=0.249;
        ZeroBoundTest=0.109;
    else
        ZeroBound=0.25/4;
        ZeroBoundTest=ZeroBound-0.01;
    end
    
    if yypred(t,ind_r)<ZeroBound && bdd_int_rate
        Shocks(t,ind_r_sh) = 0;
        Shocks(t,ind_r_sh) = (ZeroBound - DD(ind_r) - ZZ(ind_r,:)*(TTT*z + RRR*Shocks(t,:)'))/(ZZ(ind_r,:)*RRR(:,ind_r_sh));
        z = CCC+TTT*z+RRR*Shocks(t,:)';
        yypred(t,:) = (DD+ZZ*z)';
        yypred_s(t,:) = A*z;
        if yypred(t,ind_r) < ZeroBoundTest;
            error('0.25 interest rate bound procedure not working')
        end
    else
        z = z_test;
    end
    
end


