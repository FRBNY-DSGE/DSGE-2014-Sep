function [Vseq,Vseq_alt,Shseq,Shseq_alt] = setVseq(plotType,nshocks,nant,nvar,mspec,varargin)


%% Preallocate space for memory usage
Vseq=[];
Vseq_alt=[];
Shseq=[];
Shseq_alt=[];

%% Settings
switch plotType
    case {'Forecast','Shock Decomposition','Shock Decomposition-noDet','Exp_Forecast'}
      Vseq = 1:(nvar-nant-1);
        % -1 because we don't want to plot the output lev yet
      Vseq_alt = [];
    case {'Counterfactual by Variable'}
      Shseq = 1:nshocks-nant;
      Shseq_alt = [];
      Vseq = 1:(nvar-nant-1);
        % -1 because we don't want to plot the output lev yet
      Vseq_alt = [];
    case {'Counterfactual by Shock'}
      Shseq = 1:(nvar-nant-1);
        % -1 because we don't want to plot the output lev yet
      Shseq_alt = [];
      Vseq = 1:nshocks-nant; % These are now shocks
      Vseq_alt = [];
    case {'Shock','Shock Squared','Htil','Eta Squared','Sigma'}
      Vseq = 1:10; % These are now shocks
      Vseq_alt = [];
    case {'ShockAnt'}
      Vseq = [11:10+nant];
      Vseq_alt = [];
    case {'Forecast Comparison', '4Q Forecast', 'Forecast Comparison Separate', '4Q Forecast Separate'}
      Vseq = 1:nvar;
      Vseq_alt = [];
        
end

end
