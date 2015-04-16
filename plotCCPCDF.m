function [edf_res] = plot_ccp_cdf(resStruct)
%PLOT_CCP_CDF Summary of this function goes here
%   Detailed explanation goes here

% edf_res = cumsum(nanmean(resStruct.lftRes.lftHist_IIa,1)*(3/2)+nanmean(resStruct.lftRes.lftHist_Ia,1)*(3/2));
edf_res = cumsum(nanmean(resStruct.lftRes.lftHist_Ia,1)*(3));
t_res = resStruct.lftRes.t;
plot(t_res,edf_res)

end

