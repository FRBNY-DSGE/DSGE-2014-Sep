
function [nvar,varnames,graph_title,cum_for,popadj,varnames_YL,varnames_irfs,varnames_YL_4Q,varnames_YL_irfs,...
names_shocks,names_shocks_title,nl_shocks_title,shocksnames,cum_irf,vardec_varnames,shockcats,list,shockdec_color] = mspec_add(mspec,dataset,zerobound,nant,fourq_flag)

eval(['states',num2str(mspec)]);
if zerobound
    polshocks = [r_sh (r_sh+1:r_sh+nant)];
else
    polshocks = r_sh;
end

if any(mspec==[557]) 
  nvar = 7;

  varnames = {'Output Growth'; 'Aggregate Hours Growth'; 'Labor Share'; ...
              'Core PCE Inflation'; 'Interest Rate'; 'Spread'; 'Output Level'};
  varnames_irfs = varnames;

  % Save titles for graphs
  graph_title = {'Output'; 'Labor_Supply'; 'Labor_Share'; ...
                 'Inflation'; 'Interest_Rate'; 'Spread'; 'Output_Level'};
  
  % Adjustments to be made before plotting, like going to quarterly annualized
  cum_for = [1 2 3 1 0 0 2];

  % Whether to population adjust the variables
  popadj  = [1 1 0 0 0 0 1];

  % Y-Axis Labels
  varnames_YL = {'Percent Q-to-Q Annualized'; 'Percent Q-to-Q Annualized'; ... 
                  'Level'; 'Percent Q-to-Q Annualized'; 'Percent Annualized'; ...
                  'Percent Annualized'; 'Percent Annualized'};
  varnames_YL_4Q = {'Percent 4Q Growth'; 'Percent 4Q Growth'; 'Level'; ...
                    'Percent 4Q Growth'; 'Percent Annualized'; ...
                    'Percent Annualized'; 'Percent Annualized'};
  varnames_YL_irfs = {'Percent Annualized'; 'Percent Annualized'; 'Percent'; ...
                      'Percent Annualized'; 'Percent Annualized'; ...
                      'Percent Annualized'; 'Percent Annualized'};

  % Plot titles for plots of shock histories
  names_shocks = {'TFP'; 'Labor'; 'MEI'; 'b'; 'Demand'; 'Mark-Up'; 'Spread'; ...
                  '\mu_e'; '\gamma'; 'Money'};
  shocksnames = names_shocks';

  % Save names for plots of shock histories
  names_shocks_title = {'TFPShock'; 'phi'; 'mu'; 'b'; 'g'; 'lambda_f'; 'w'; ...
                        'mu_e'; 'gamma'; 'Money'};
                                         
             
  nl_shocks_title = {'Technology'; '\phi'; 'Financial'; 'b'; 'g'; 'Mark-Up'; 'w'; ...
                     '\mu_e'; '\gamma'; 'Monetary Policy'};
                                         
  % Not relevant
  cum_irf = [];
  vardec_varnames = {'Output Growth'; 'Aggregate Hours Growth'; 'Labor Share';...
                    'Core PCE Inflation'; 'Interest Rate'; 'Spread'; };

  % Which shocks to plot in the shock decompositions: 
  % - Give the indices of the shocks, where the labels (used below) aliasing
  %   the index numbers come from running statesMSPEC above
  % - Can group shocks by making an entry be an array like [mu_sh z_sh]
  shockcats = {sigw_sh; mu_sh; z_sh; polshocks; laf_sh; g_sh; phi_sh;};
    % shockdogs -- deprecated because it didn't work well with shockcats

  % How to label each element of shockcats in the legend
  list = {'Spread', 'MEI', 'TFP', 'Policy', 'Mark-Up', 'Gov''t', 'Labor'};

  % Colors of the bars in the shockdecs
  shockdec_color = {'indigo','aqua','firebrick','darkorange','yellowgreen','gold','pink'};

  %% Adding information for expectatons
  if zerobound 
    gen_temp = @(pattern) arrayfun(@(i_) sprintf(pattern, i_), 1:nant, 'UniformOutput', false);
    temp1 = gen_temp('exp_%i');
    temp2 = gen_temp('');
    temp3 = gen_temp('Ant %i');
    temp4 = gen_temp('%iqA ant.sh.');
    temp5 = gen_temp('Ant%i');

    varnames(end+1:end+nant)           = temp1(1:nant);
    varnames_YL(end+1:end+nant)        = temp2(1:nant);
    varnames_YL_4Q(end+1:end+nant)     = temp2(1:nant);
    varnames_irfs(end+1:end+nant)      = temp2(1:nant);
    varnames_YL_irfs(end+1:end+nant)   = temp2(1:nant);
    cum_for(end+1:end+nant)            = zeros(1,nant);
    popadj(end+1:end+nant)             = zeros(1,nant);
    names_shocks(end+1:end+nant)       = temp3(1:nant);
    nl_shocks_title(end+1:end+nant)    = temp4(1:nant);
    names_shocks_title(end+1:end+nant) = temp5(1:nant);

    shocksnames = names_shocks';
  end

end

