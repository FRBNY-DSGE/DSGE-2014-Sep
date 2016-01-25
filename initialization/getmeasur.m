%%
% getmeasur.m: This is a driver for the mspec-specific measur function.
%
function [ZZ,DD,DDcointadd,QQ,EE,MM,retcode] = ...
  getmeasur(mspec,TTT,RRR,valid,params,nvar,nlags,npara,coint,cointadd,varargin)

  if nargin > 10
    extra_arg = varargin(1);
  else
    extra_arg = {};
  end

  [ZZ,DD,DDcointadd,QQ,EE,MM,retcode] = feval(['measur' num2str(mspec)], ...
    TTT,RRR,valid,params,nvar,nlags,mspec,npara,coint,cointadd,extra_arg{:});

end
