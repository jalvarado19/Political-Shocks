% Computes Implulse Response Functions to a Monetary Policy Shocks
% (identified with Cholesky orthogonalization, ordering the FFR last)
% using the BVAR of Giannone, Lenza and Primiceri (2012)
% Based on the default setting, including all priors
% %%%%%%%%%%%%%%%%%%%%
% The model includes the following data 
% RGDP: 4 x logarithm of Real Gross Domestic Product, Quantity Index (2000=100) , SAAR
% PGDP: 4 x logarithm of Gross domestic product Price Index
% Cons: 4 x logarithm of Real Personal Consumption Expenditures, Quantity Index (2000=100) , SAAR
% GPDInv: 4 x logarithm of Real Gross Private Domestic Investment, Quantity Index (2000=100) , SAAR
% Emp. Hours: 4 x logarithm of HOURS OF ALL PERSONS: NONFARM BUSINESS SEC (1982=100,SA)
% Real Comp/Hour: 4 x logarithm of REAL COMPENSATION PER HOUR,EMPLOYEES:NONFARM BUSINESS(82=100,SA)
% FedFunds: INTEREST RATE: FEDERAL FUNDS (EFFECTIVE) (% PER ANNUM,NSA)
% 
% clear all
% close all
% addpath([cd '/subroutines'])  %on a MAC
% %addpath([cd '\subroutines']) %on a PC
% 
% load DataSW 
% Load data from the dataset of Stock and Watson (2008)
% The variables enter the models in annualized log-levels (i.e. we take logs and multiply by 4), 
% except those already defined in terms of annualized rates, such as
% interest rates, which are taken in levels and are divided by 100.



lags = 5;

% Run the Bayesian VAR
res = bvarGLP(Y1,lags,'mcmc',1,'pos',8,'MCMCconst',1);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Arrange and display some results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Computes the Impulse response function

hmax =100; % maximum horizon for the IRFs
nshock = 8; % Position of the variable to shock; 7 is the FFR

% IRfs at the posterior mode 
beta = res.postmax.betahat;
sigma = res.postmax.sigmahat;
irf =  bvarIrfs(beta,sigma,nshock,hmax);


% IRFs at each draw
ndraws = size(res.mcmc.beta,3);

Dirf = zeros(hmax,size(Y1,2),ndraws);

for jg = 1:ndraws
    beta  = res.mcmc.beta(:,:,jg);
    sigma = res.mcmc.sigma(:,:,jg);
    Dirf(:,:,jg) =  bvarIrfs(beta,sigma,nshock,hmax);
end

sIRF = sort(Dirf,3);


%plots the IRFs to a Monetary Policy Shock
ShortDescr{1} = 'Real GDP';
ShortDescr{2} = 'I';
ShortDescr{3} = 'TFP';
ShortDescr{4} = 'G';
ShortDescr{5} = 'Tax';
ShortDescr{6} = 'Transf';
ShortDescr{7} = 'Sp500';
ShortDescr{8} = 'Ideology';

%%
for jn = 1:8
    
    if jn <8
        
        subplot(2,4,jn)
        plot(0:hmax-1,irf(:,jn)/4*100,0:hmax-1,squeeze(sIRF(:,jn,round([.16 .5 .84]*ndraws)))/4*100,'-.r')
        title(ShortDescr{jn})
    
    else
        
        subplot(2,4,jn)
        plot(0:hmax-1,irf(:,jn)*100,0:hmax-1,squeeze(sIRF(:,jn,round([.16 .5 .84]*ndraws)))*100,'-.r')
        title(ShortDescr{jn})
        
        legend('IRF at mode','16th quantile', '50th  quantile', '84th quantile','Location','NorthEastOutside')
    
    end
    
end

