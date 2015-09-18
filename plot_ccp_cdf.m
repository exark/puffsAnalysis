function [edf_res, t_res] = plot_ccp_cdf(resStruct,color)
% plot_ccp_cdf prints a graph of the cumulative lifetime distribution of
% the passed resStruct, in the desired color.

summedHist = sum(resStruct.lftRes.lftHist_total);
summedHist = summedHist/sum(summedHist);
edf_res = cumsum(summedHist);
t_res = resStruct.lftRes.t;
plot(t_res,edf_res,'color',color);

end