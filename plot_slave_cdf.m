function [edf_res, t_res] = plot_slave_cdf(resStruct,color)
% plot_slave_cdf prints a graph of the cumulative lifetime distribution of
% CCPs which contain the slave from resStruct, in the desired color.

slaveHist = resStruct.lftRes.lftHistSlaveAll{1,1};
edf_res = nanmean(cumsum(slaveHist*3,2));
t_res = resStruct.lftRes.t;
plot(t_res,edf_res,'color',color);

end