function [zprev, pprev] = augmentStates(mspec,nant,nstate,zend,Pend)

% Because money shock named differently in different models
if mspec == 904
    r_tl1 = getState(mspec,nant,'rm_tl1');
else
    r_tl1 = getState(mspec,nant,'r_tl1');
end
r_tlx = r_tl1;
zprev = zeros(nstate,1);
zprev(1:r_tl1-2) = zend(1:r_tlx-2,end);
zprev(r_tl1+nant:end) = zend(r_tlx-1:end,end);

pprev = zeros([nstate nstate]);
pprev(1:r_tlx-2,1:r_tlx-2) = Pend(1:r_tlx-2,1:r_tlx-2,end);
pprev(r_tl1+nant:end,1:r_tlx-2) = Pend(r_tlx-1:end,1:r_tlx-2,end);
pprev(1:r_tlx-2,r_tl1+nant:end) = Pend(1:r_tlx-2,r_tlx-1:end,end);
pprev(r_tl1+nant:end,r_tl1+nant:end) = Pend(r_tlx-1:end,r_tlx-1:end,end);