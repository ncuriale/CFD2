function plotSol( P_d, M_d, rho_d,P_exact, M_exact, RHO_exact, bcFlag, nodes)

if (bcFlag == 1)
  
    %Mach
    figure(1)
    axes('YMinorTick','on','YGrid','on','XMinorTick','on','XGrid','on');
    xlim([0 nodes]);
    box('on');
    hold('all');
    plot(M_exact(2:nodes),'r-','LineWidth',1.5);
    plot(M_d(1:nodes-1),'kx');
    grid on
    xlabel('Nodes','interpreter','latex','FontSize',13);
    ylabel('Mach Number (-)','interpreter','latex','FontSize',13);
    legend({'Exact','Non-Diagonal'},'interpreter','latex','FontSize',11,'Location','NorthEast')
    title('Subsonic Nozzle -- 99 nodes, Cn=120, k=(0.5,0.02)','interpreter','latex','FontSize',12);
       
    %Pressure
    figure(2)
    axes('YMinorTick','on','YGrid','on','XMinorTick','on','XGrid','on');
    xlim([0 nodes]);
    box('on');
    hold('all');
    plot(P_exact(2:nodes),'r-','LineWidth',1.5);
    plot(P_d(1:nodes-1),'kx');
    grid on
    xlabel('Nodes','interpreter','latex','FontSize',13);
    ylabel('Pressure (Pa)','interpreter','latex','FontSize',13);
    legend({'Exact','Non-Diagonal'},'interpreter','latex','FontSize',11,'Location','NorthEast')
    title('Subsonic Nozzle -- 99 nodes, Cn=120, k=(0.5,0.02)','interpreter','latex','FontSize',12);

elseif (bcFlag == 2) 

    %Mach
    figure(1)
    axes('YMinorTick','on','YGrid','on','XMinorTick','on','XGrid','on');
    xlim([0 nodes]);
    box('on');
    hold('all');
    plot(M_exact(2:nodes),'r-','LineWidth',1.5);
    plot(M_d(1:nodes-1),'kx');
    grid on
    xlabel('Nodes','interpreter','latex','FontSize',13);
    ylabel('Mach Number (-)','interpreter','latex','FontSize',13);
    legend({'Exact','Non-Diagonal'},'interpreter','latex','FontSize',11,'Location','NorthWest')
    title('Transonic Nozzle -- 99 nodes, Cn=40, k=(0.0,0.02)','interpreter','latex','FontSize',12);
        
    %Pressure
    figure(2)
    axes('YMinorTick','on','YGrid','on','XMinorTick','on','XGrid','on');
    xlim([0 nodes]);
    box('on');
    hold('all');
    plot(P_exact(2:nodes),'r-','LineWidth',1.5);
    plot(P_d(1:nodes-1),'kx');
    grid on
    xlabel('Nodes','interpreter','latex','FontSize',13);
    ylabel('Pressure (Pa)','interpreter','latex','FontSize',13);
    legend({'Exact','Non-Diagonal'},'interpreter','latex','FontSize',11,'Location','NorthEast')
    title('Transonic Nozzle -- 99 nodes, Cn=40, k=(0.0,0.02)','interpreter','latex','FontSize',12);
  
elseif (bcFlag == 3)

    %Mach
    figure(1)
    axes('YMinorTick','on','YGrid','on','XMinorTick','on','XGrid','on');
    xlim([0 nodes]);
    box('on');
    hold('all');
    plot(M_exact(2:nodes),'r-','LineWidth',1.5);
    plot(M_d(1:nodes-1),'k-x');
    grid on
    xlabel('Nodes','interpreter','latex','FontSize',13);
    ylabel('Mach Number (-)','interpreter','latex','FontSize',13);
    legend({'Exact','Non-Diagonal'},'interpreter','latex','FontSize',11,'Location','NorthEast')
    title('Shock Tube @ t=0.0061ms -- 399 nodes, Cn=1, k=(0.5,0.02)','interpreter','latex','FontSize',12);
        
    %Density
    figure(2)
    axes('YMinorTick','on','YGrid','on','XMinorTick','on','XGrid','on');
    xlim([0 nodes]);
    box('on');
    hold('all');
    plot(RHO_exact(2:nodes),'r-','LineWidth',1.5);
    plot(rho_d(1:nodes-1),'k-x');
    grid on
    xlabel('Nodes','interpreter','latex','FontSize',13);
    ylabel('Density (kg/m3)','interpreter','latex','FontSize',13);
    legend({'Exact','Non-Diagonal'},'interpreter','latex','FontSize',11,'Location','NorthEast')
    title('Shock Tube @ t=0.0061ms -- 399 nodes, Cn=1, k=(0.5,0.02)','interpreter','latex','FontSize',12);
      
end
    

