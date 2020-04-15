%% How much politics matter for Macro: Evidence from US

%% Data
clear
close all

run data_macro

%% Save date

%Y1 = [Y_47(:,1) Y_47(:,2) Y_47(:,3)  Y_47(:,5) Y_47(:,6) Y_47(:,4) Y_47(:,7)  [Pol(2:end,22);Pol(end,22)]];
%Y1 = [Y_47(:,1) Y_47(:,2) Y_47(:,3)  Y_47(:,5) Y_47(:,6) Y_47(:,4) Y_47(:,7)  Pol(1:end,22)];

Y1  = [Y_1960(1:end-1,1) Y_1960(1:end-1,2) Y_1960(1:end-1,3) Y_1960(1:end-1,4) Y_1960(1:end-1,5) Y_1960(1:end-1,6) log(Y_1960(1:end-1,11))  Pol(53:end,22)];
%%
tic

%% Model 1 Reduced form
%%

model = 3;

run IRFs 

if model == 3
    
global n p T
  
  [T,n]=size(Y1);
  p = lags;
  show_graph = 0;
  res = bvarGLP(Y1,lags,'mcmc',0,'MCMCconst',1);

  x0     = zeros(57,1);
  
  x0(1)  = 1;
  x0(10) = 1;
  x0(19) = 1;
  x0(28) = 1;
  x0(37) = 1;
  x0(46) = 1;
  x0(55) = 1;
  x0(57) = 1;
  
  beta_aux = res.postmax.betahat;
  b_prior  = res.postmax.bprior;
  S        = (res.postmax.epshat)'*(res.postmax.epshat);
  sigmahat = res.postmax.sigmahat;
  Omega    = diag(res.postmax.omega);
  V_aux    = res.postmax.V;
  
  [fh,A,gh,H,itct,fcount,retcodeh] = csminwel(@(x) posterior_A(x,beta_aux,b_prior,S,Omega),x0(:),eye(57)*1.1,'',10e-5,1000);
  
    A0 = construct_A(A);
    
    Gamma = inv(A0);
    
    hmax =8; % maximum horizon for the IRFs
    nshock = 8;
    
    %% MCMS Algorithm
    
    it = 2;
    maxit = 20000;
    
    a0 = A;
    
    a_test1 = zeros(1,maxit);
    a_test2 = zeros(1,maxit);
    a_test3 = zeros(1,maxit);
    
    a_test1(1)=a0(1);
    a_test2(1)=a0(12);
    a_test3(1)=a0(18);
    
    A0_draw(:,1) = a0;
    
    c=0.5;
    nreject  = 0;
    nrejec2  = 0;
    
    while it<maxit
        it
        R = mvnrnd(a0,c^2*H);
        
         A0_aux = construct_A(R);

        alpha = min(1,abs(exp(posterior_A2(R,beta_aux,b_prior,S,Omega)-posterior_A2(a0',beta_aux,b_prior,S,Omega))));
        reject = rand(1);
        

        if reject < alpha
            struct_shock = A0_aux*(res.postmax.epshat)';
            
            struct_shock = struct_shock(:,1:T-p)';
            
            elections = repmat(1-election_date(58:end),1,n);
            
            impact = (abs(struct_shock).*elections);
            
            nrejec2  = nrejec2 +1;
            %shock = sum(impact(:,8));
            shock = sum(impact(:,8));
            %0.29485
            if shock/T<0.289 && A0_aux(n,n)>0
                a0 = R';
                a_test1(it) = R(1);
                a_test2(it) = R(12);
                a_test3(it) = R(18); 
                nreject = nreject+1;
                A0_draw(:,it) = R';
            end
            
        else
            A0_draw(:,it) = a0;
            a_test1(it) = a0(1);
            a_test2(it) = a0(12);
            a_test3(it) = a0(18);
        end 
        it=it+1;
    end
    rate1 = nreject/maxit
    rate2 = nrejec2/maxit
%%
 run graphs_mcs
 

    hmax =16; % maximum horizon for the IRFs
    nshock = 8; % Position of the variable to shock; 7 is the FFR

    % IRfs at the posterior mode 
    beta = res.postmax.betahat;
    sigma = res.postmax.sigmahat;
    irf =  bvarIrfs(beta,sigma,nshock,hmax,A0,1);


    % IRFs at each draw
    ndraws = maxit; %size(res.mcmc.beta,3);
    burn = 10000;
    Dirf = zeros(hmax,size(Y1,2),ndraws-burn);

    for jg = 1:ndraws-burn-1
%         beta  = res.mcmc.beta(:,:,jg);
%         sigma = res.mcmc.sigma(:,:,jg);
         Bsample = mvnrnd(beta_aux(:),kron((inv(construct_A(A0_draw(:,burn-1+jg)))*inv(construct_A(A0_draw(:,burn-1+jg)))'),V_aux));
         B_mat   = reshape(Bsample,[41,8]);
        
        Dirf(:,:,jg) =  bvarIrfs(B_mat,sigma,nshock,hmax,construct_A(A0_draw(:,burn-1+jg)),1);
    end

    sIRF = sort(Dirf,3);
    
    test = squeeze(sIRF(:,8,:));


    %plots the IRFs to a Monetary Policy Shock
    ShortDescr{1} = 'Real GDP';
    ShortDescr{2} = 'Investment';
    ShortDescr{3} = 'TFP';
    ShortDescr{4} = 'G';
    ShortDescr{5} = 'Tax revenue';
    ShortDescr{6} = 'Transfers';
    ShortDescr{7} = 'S&P500';
    ShortDescr{8} = 'Political Shock';

  
    for jn = 1:8

        if jn <8

            subplot(2,4,jn)
            plot(0:hmax-1,irf(:,jn)*100,0:hmax-1,squeeze(sIRF(:,jn,round([0.15 0.5 .85]*(ndraws-burn)))*100),'LineWidth',2)
            title(ShortDescr{jn})

        else

            subplot(2,4,jn)
            plot(0:hmax-1,irf(:,jn)*100,0:hmax-1,squeeze(sIRF(:,jn,round([.15 .5 .85]*(ndraws-burn)))*100),'LineWidth',2)
            title(ShortDescr{jn})


        end

    end
end

%%
if model == 4
    
  res = bvarGLP(Y1,lags,'mcmc',0,'MCMCconst',1);
  x0 = zeros(22,1);
  
  x0(1)=1;
  x0(3)=1;
  x0(6)=1;
  x0(7)=1;
  x0(10)=1;
  x0(13)=1;
  x0(20)=1;
  x0(22)=1;
  
  
  
  beta_aux = res.postmax.betahat;
  b_prior  = res.postmax.bprior;
  S        = (res.postmax.epshat)'*(res.postmax.epshat);
  sigmahat = res.postmax.sigmahat;
  Omega    = diag(res.postmax.omega);
  
%   A0 = [A(1) 0     0    0     0     0     0    0;
%         A(2) A(3)  0    0     0     0     0    0;
%         A(4) A(5)  A(6) 0     0     0     0    0;
%         0    0     0    A(7)  A(8)  0     0    0;
%         0    0     0    A(9)  A(10) 0     0    0;
%         0    0     0    A(11) A(12) A(13) 0    0;
%         A(14) A(15)  A(16) A(17) A(18) A(19)  A(20) A(21);
%         0    0     0    0    0     0     0    A(22)];

  
  global n p T
  
  [T,n]=size(Y1);
  p = lags;
  
  [fh,A,gh,H,itct,fcount,retcodeh] = csminwel(@(x) posterior_A(x,beta_aux,b_prior,S,Omega),x0(:),eye(22)*1.01,'',10e-5,1000);
  
    A0 = [A(1) 0     0    0     0     0     0    0;
        A(2) A(3)  0    0     0     0     0    0;
        A(4) A(5)  A(6) 0     0     0     0    0;
        0    0     0    A(7)  A(8)  0     0    0;
        0    0     0    A(9)  A(10) 0     0    0;
        0    0     0    A(11) A(12) A(13) 0    0;
        A(14) A(15)  A(16) A(17) A(18) A(19)  A(20) A(21);
        0    0     0    0    0     0     0    A(22)];
    
    hmax =8; % maximum horizon for the IRFs
    nshock = 8;
    %%
    it = 2;
    maxit = 20000;
    
    a0 = A;
    
    a_test1=zeros(1,maxit);
    a_test2=zeros(1,maxit);
    a_test3=zeros(1,maxit);
    
    a_test1(1)=a0(1);
    a_test2(1)=a0(12);
    a_test3(1)=a0(18);
    
    a_test12=zeros(1,maxit);
    a_test22=zeros(1,maxit);
    a_test32=zeros(1,maxit);

    
    A0_draw(:,1) = a0;
    
    c=40;
    nreject=0;
    nreject2=0;
    while it<maxit
        R = mvnrnd(a0,c^2*H);
        R2 = mvnrnd(a02,c^2*H);
        
        alpha = min(1,abs(posterior_A2(R,beta_aux,b_prior,S,Omega)/posterior_A2(a0,beta_aux,b_prior,S,Omega)));
        reject = rand(1);
        
        alpha2 = min(1,abs(posterior_A2(R2,beta_aux,b_prior,S,Omega)/posterior_A2(a02,beta_aux,b_prior,S,Omega)));
        reject2 = rand(1);
        
        if reject < alpha
            a0 = R';
            a_test1(it) = R(1);
            a_test2(it) = R(12);
            a_test3(it) = R(18);
            A0_draw(:,it) = R';
            nreject = nreject+1;
        else
            a_test1(it) = a0(1);
            a_test2(it) = a0(12);
            a_test3(it) = a0(18);
            A0_draw(:,it) = a0;
        end 
        
        it=it+1;

    end
    
    %%
    
    irf =  bvarIrfs(beta_aux,sigmahat,nshock,hmax,A0,1);
     
end


