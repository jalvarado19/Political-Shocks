function [res1] = acceptance(c,A,beta_aux,b_prior,S,Omega)
    
    it = 2;
    maxit = 200000;
    
    a0 = A;
    
    a_test1 = zeros(1,maxit);
    a_test2 = zeros(1,maxit);
    a_test3 = zeros(1,maxit);
    
    a_test1(1)=a0(1);
    a_test2(1)=a0(12);
    a_test3(1)=a0(18);
    
    A0_draw = zeros(n*n,maxit);
    A0_draw(:,1) = a0;
    

    nreject  = 0;
    nrejec2  = 0;
    
    while it<maxit
        it
        R = mvnrnd(a0,c^2*H);
        
        A0_aux = vec2mat(R,n);
        
        alpha = min(1,abs(posterior_A2(R ,beta_aux,b_prior,S,Omega)/posterior_A2(a0',beta_aux,b_prior,S,Omega,tol)));
        reject = rand(1);
        

        if reject < alpha
            struct_shock = A0_aux*(res.postmax.epshat)';
            
            struct_shock = struct_shock(:,1:T-p)';
            
            elections = repmat(1-election_date(6:end),1,n);
            
            impact = (abs(struct_shock).*elections);
            
            nrejec2  = nrejec2 +1;
            if sum(impact(:,8))/T<tol
                a0 = R';
                a_test1(it) = R(1);
                a_test2(it) = R(12);
                a_test3(it) = R(18);
                A0_draw(:,it) = R';
                nreject = nreject+1;
            end
            
        else
            a_test1(it) = a0(1);
            a_test2(it) = a0(12);
            a_test3(it) = a0(18);
            A0_draw(:,it) = a0;
        end 
        it=it+1;
    end
    res1 = (nreject/maxit-0.25)^2;
    res2 = (nrejec2/maxit-0.25)^2;
end

