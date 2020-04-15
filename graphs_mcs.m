if show_graph ==1
   Lagmax = 1000;
   rho_n1=zeros(1,Lagmax);
   rho_n2=zeros(1,Lagmax);
   rho_n3=zeros(1,Lagmax);
    for h=1:Lagmax
        p1 = corrcoef([a_test1(h+1:end)' a_test1(1:end-h)']);
        rho_n1(h)= p1(1,2);
        p2 = corrcoef([a_test2(h+1:end)' a_test2(1:end-h)']);
        rho_n2(h)= p2(1,2);
        p3 = corrcoef([a_test3(h+1:end)' a_test3(1:end-h)']);
        rho_n3(h)= p3(1,2);
    end
    
    
    acceptance_rate = nreject/maxit*100
    
    
    figure;
    subplot(3,1,1)
    hold on
    plot(a_test1)
    hold off
    subplot(3,1,2)
    hold on
    plot(a_test2)
    hold off
    subplot(3,1,3)
    hold on
    plot(a_test3)
    hold off
    
        
    figure;
    subplot(3,6,1)
    plot(A0_draw(1,:))
    subplot(3,6,2)
    plot(A0_draw(2,:))
    subplot(3,6,3)
    plot(A0_draw(3,:))
    subplot(3,6,4)
    plot(A0_draw(4,:))
    subplot(3,6,5)
    plot(A0_draw(5,:))
    subplot(3,6,6)
    plot(A0_draw(6,:))
    subplot(3,6,7)
    plot(A0_draw(7,:))
    subplot(3,6,8)
    plot(A0_draw(8,:))
    subplot(3,6,9)
    plot(A0_draw(9,:))
    subplot(3,6,10)
    plot(A0_draw(10,:))
    subplot(3,6,11)
    plot(A0_draw(11,:))
    subplot(3,6,12)
    plot(A0_draw(12,:))
    subplot(3,6,13)
    plot(A0_draw(13,:))
    subplot(3,6,14)
    plot(A0_draw(14,:))
    subplot(3,6,15)
    plot(A0_draw(15,:))
    subplot(3,6,16)
    plot(A0_draw(16,:))
    subplot(3,6,17)
    plot(A0_draw(17,:))
    subplot(3,6,18)
    plot(A0_draw(18,:))
   
        
    figure;
    subplot(3,1,1)
    plot(rho_n1)
    subplot(3,1,2)
    plot(rho_n2)
    subplot(3,1,3)
    plot(rho_n3)
    hold off
 end