hold on
plot(Pol(:,1)*0.1)
plot(Pol(:,16))
hold off




%% Congress Statistics
mean_ideo_d = sum((Pol(:,1).*Pol(:,22)))/sum((Pol(:,1)));
mean_ideo_r = sum((1-Pol(:,1)).*Pol(:,22))/sum(1-Pol(:,1));
mean_ideo_d_full = sum((Pol(:,4).*Pol(:,1).*Pol(:,14)))/sum((Pol(:,1)).*Pol(:,4));
mean_ideo_r_full = sum((Pol(:,4).*(1-Pol(:,1)).*Pol(:,14)))/sum((1-Pol(:,1)).*Pol(:,4));

g_gdp = exp(Y_47(5:end,1))./exp(Y_47(1:end-4,1))-1;

mean_growth_d = sum((Pol(5:end,1).*g_gdp))/sum((Pol(5:end,1)));
mean_growth_r = sum((1-Pol(5:end,1)).*g_gdp)/sum((1-Pol(5:end,1)));

mean_growth_d_full = sum((Pol(5:end,4).*Pol(5:end,1).*g_gdp))/sum((Pol(5:end,1)).*Pol(5:end,4));
mean_growth_r_full = sum((Pol(5:end,4).*(1-Pol(5:end,1)).*g_gdp))/sum((1-Pol(5:end,1)).*Pol(5:end,4));

mean(g_gdp)

hold on
plot(g_gdp,'LineWidth',2)
plot(Pol(5:end,22)/10,'LineWidth',2)
hold off

hold on
plot(Y_1960(:,11)/std(Y_1960(:,11)))
plot(Pol(52:end,22)/std(Pol(52:end,22)))
hold off

corrcoef(g_gdp,Pol(5:end,22))
corrcoef(Pol(:,1),Pol(:,16))
corrcoef(Y_1960(:,11),Pol(52:end,22))

ols1 = fitlm([Pol(5:end,1) Pol(5:end,22)], g_gdp,'Intercept',false);
ols2 = fitlm(Pol(5:end,22), g_gdp,'Intercept',false);
ols3 = fitlm([Pol(5:end,1) Pol(5:end,16) (Pol(5:end,4).*Pol(5:end,1).*Pol(5:end,16)) Pol(5:end,4).*(1-Pol(5:end,1)).*Pol(5:end,16)], g_gdp,'Intercept',false);
ols4 = fitlm([Pol(5:end,1) Pol(5:end,16) (Pol(5:end,4).*Pol(5:end,1).*Pol(5:end,16)) Pol(5:end,4).*(1-Pol(5:end,1)).*Pol(5:end,16)], g_gdp,'Intercept',true);


