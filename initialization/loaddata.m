%% loaddata
%
% Description: reads time series data from ASCII 
% Output: 
%        1) YYall ---> matrix of observables
%        2) XXall ---> 
%        3) ti --> 
%        4) nobs --> number of periods in data imported
%        4) dlpopall --> log differences of the population obtained from Haver Analaytics
%        5) dlMA_pop --> log differences of the population forecasted by Macro Advisers 

function [YYall,XXall,ti,nobs,dlpopall,dlMA_pop,MA_pop,population] = loaddata(nvar,nlags,nant,antlags,psize,zerobound,peachflag)

    
start_date = 1959.00;
load('rawData');
    
    popreal=data(:,10);
    m = length(MA_pop);
    [xxx,popfor]=Hpfilter([popreal;MA_pop(2:end)],1600);
    population = popfor(1:end-m+1);
    MA_pop = popfor(end-m+1:end);
    
    dlpopreal = log(popreal(2:end)) - log(popreal(1:end-1));
    dlpopall = log(population(2:end)) - log(population(1:end-1));
    dlMA_pop = log(MA_pop(2:end)) - log(MA_pop(1:end-1));
    
    
    nobs   = size(data,1);
    ti     = (start_date:0.25:(start_date+.25*(nobs-1)))';
    
    %% Compute growth rates
    series  = zeros(nobs-1,nvar);
    nobs = nobs - 1;
    ti = ti(2:end);

    %% Output growth (log approximation quarterly annualized
    %% Load levels
    Output =  400*(log(data(2:end,1))-log(data(1:end-1,1)));    
    series(:,1) = Output + 400*(dlpopreal - dlpopall);
    %% lhoursuppc: log hours worked per capita
     hours = 100*log(data(2:end,4));

     series(:,2) = hours + 100*(log(popreal(2:end))-log(population(2:end)));
    
    %% log labor income share
    series(:,3) = 100*log(data(2:end,6));
    
    %% Change to core PCE - annualized
    series(:,4) = 400*(log(data(2:end,7))-log(data(1:end-1,7)));  
    
    %% nominal short-term interest rate (3 months) - % annualized
    series(:,5) = data(2:end,8);
    
    %% spread: BAA-10yr TBill
    series(:,6) = data(2:end,12);

    series(:,7) = 100*log(data(2:end,1)) + 100*(log(popreal(2:end))-log(population(2:end)));

    [ExpFFR, peachdata_FFR] = ExpFFR_OIS(nant,antlags,psize,zerobound,peachflag);
    E_num = size(ExpFFR,1);
    fill = nan(length(series)-E_num,nant);
    data_add = [fill;ExpFFR];
    series = [series,data_add];
    nvar = nvar + nant;
    
    YYall     = series(1+nlags:nobs,:);%%nlags switched from T0, DF
    XXall     = ones(nobs-nlags, 1+nvar*nlags);
    ti     = ti(1+nlags:nobs,:);
    dlpopall = dlpopall(1+nlags:nobs); % NEW
    
    
    for i = 1:1:nlags;
      XXall(:,1+(i-1)*nvar+1:1+i*nvar) = series(nlags-i+1:nobs-i,:);
    end
    
    if nlags > 0
        XXall = [XXall;[1,YYall(end,:),XXall(end,1+1:1+nvar*nlags-nvar)]];
    end
    nobs = size(YYall,1);


%% Plot Input Data

    datt = [ti,YYall];
    FORMAT = '\n %2.2f  ';
    
    for i = 1:nvar
        FORMAT = strcat(FORMAT,' %2.2f ');
    end
    fprintf(1,'\n\n Actual data  \n\n');
    fprintf(1,FORMAT,datt');
    fprintf('\n')
    
    for I = 1:size(YYall,2)
        plot(ti,YYall(:,I))
        close all;
    end

