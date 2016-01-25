% OVERVIEW
% dsgesolv.m finds a solution to the DSGE model with an input parameter vector.
%            To help construct the G0, G1, C, PSI, and PIE matrices, which are used to
%            express the model in canonical form, the following programs are called
%            (1) getpara<>.m: assigns a parameter name to each element of a
%                            parameter vector. Also assigns values to steady state parameters
%                            (parameters that are derived from the values of other parameters).
%                            These paramters are usually used in the steady state component
%                            (long-term trend, D) of the measurement equation: Y_t = Z*S_t +D
%            (2) states<mspec>.m: Assigns a number to each state.
%            (3) eqs<mspec>.m: Assigns a number to each equation.
%            (4) eqcond<mspec>.m: Expresses the equilibrium conditions in
%                                canonical form using G0, G1, C, PSI, and PIE matrices.
%
%                                Using the assigned states and equations in states*.m and eqs*.m,
%                                coefficients are specified in the proper positions.
%                                G0 (neq x neq) holds coefficients of current time states.
%                                G1 (neq x neq) holds coefficients of lagged states.
%                                PSI (neq x nshocks) holds coefficients of the state corresponding to an iid shock.
%                                PIE (neq x n expecational states) holds coefficients of expectational states.
% 
% INPUTS
% para: a parameter vector
% mspec: model specification number
% nant: number of anticipated shocks
% varargin: Optional non-zero integer value used to specify 
%           alternative policy rules, which are 
%           written in eqcond file of the model
%   ** 0 should be reserved for "No Changes/Regular Rule"

% OUTPUTS
% TTT,RRR,CCC: matrices of tht state transition equation:
%              S_t = TTT*S_(t-1) + RRR*eps_t + CCC
% valid: If the solution to the model meets any of the various error
%        critera in this program, this flag indicates the type of error
%        encountered. 
%
% See also: Sims(2000), Solving Linear Rational Expectations Models

function [TTT,RRR,CCC,valid] = dsgesolv(mspec, para, nant)

valid = 1;
TTT = 0;
RRR = 0;
CCC = 0;

%% Assign a value to each parameter
getPara_script

%% Pre-allocate matrices


if exist('fzflag','var')
    if fzflag < 1
        disp('fzero failed to converge');
        TTT = [];
        RRR = [];
        valid = -1;
        return
    end
end

%% Assign a number to each state

eval(['states',num2str(mspec)]);

%% Assign numbers to each equations

eval(['eqs',num2str(mspec)]);

%% Express the equilibrium conditions in canonical form

% Total number of equations, states, shocks, and endogenous variables

    nstate = n_end+n_exo+n_exp;
    neq  = nstate;

G0 = zeros(neq,neq);
G1 = zeros(neq,neq);
C =  zeros(neq,1);
PSI = zeros(neq,nex);
PIE =  zeros(neq,nend);


eval(['eqcond',num2str(mspec)]);

if any(isnan(G0(:))) || any(isnan(G1(:)))
    disp('NaN values in G0 or G1');
    %keyboard;
    valid = -2;
    return
end


try
    [T1,TC,T0,M,TZ,TY,gev,RC] = gensys(G0,G1,C,PSI,PIE,1+10^(-6));
catch err
    warning('gensys failed to run properly');
    T1 = zeros(nstates,nstates);
    T0 = zeros(nstates,2);
    TC = zeros(nstates,1);
    RC = [1,1];
    valid = -4;
end

if (RC(1) ~= 1) || (RC(2) ~= 1);
    TTT = zeros(nstates,nstates);
    RRR = zeros(nstates,2);
    valid = - 3;
    %keyboard;
end;

TTT =    real( T1 );
RRR =    real( T0 );
CCC =    TC;

% Some of our observables are growth rates, which is calculated as a 
% linear combination of a present and lagged state. To capture the lagged state, 
% we assign to it an index. In addition, we also need to expand the 
% matrices of the state transition equation to accommodate the extra state. 
% In dsgesolv.m, AddTTT is appended to TTT to capture the lagged state 
% value in the current state vector, while AddRRR and AddCCC augment 
% RRR and CCC with the appropriate number of zeros.
    numAdd = 1;
    AddTTT = zeros(numAdd,size(TTT,2));
    AddTTT(1,y_t) = 1;
    TTT = [[TTT,zeros(size(TTT,1),numAdd)];[AddTTT,zeros(numAdd)]];
    RRR = [RRR ; zeros(numAdd,size(RRR,2))];
    CCC = [CCC ; zeros(numAdd,size(CCC,2))];

    

