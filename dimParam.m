function [P_d, rho_d, u_d, M_d, T_d, c_d, e_d] = dimParam(S, Q, nodes, bcFlag)

% constants
R = 287; 
gam = 1.4; 

% Boundary conditions
if (bcFlag == 1) 

    uL = 65.4510; 
    PL = 97534.31;
    rhoL = 1.1409;
    eL = 215864.21;%PLo/(gam-1) + 0.5*rhoLo*uLo^2;%2443.63 215864.21
    cL = sqrt(gam*PL/rhoL); 
    SL = 2.5; 
    QL =[rhoL*SL; rhoL*uL*SL; rhoL*eL*SL];
    uR = 113.0563;
    PR = 92772.07; 
    rhoR = 1.1008; 
    eR = 217083.84;%PRo/(gam-1) + 0.5*uRo^2; 
    cR = sqrt(gam*PR/rhoR);
    SR = 1.5;
    QR =[rhoR*SR; rhoR*uR*SR; rhoR*eR*SR]; 

%     % Non-dimensionalized 
%     uL = uLo/cLo; 
%     rhoL = rhoLo/rhoLo;
%     eL = (rhoLo*ELo)/(rhoLo*cLo^2);
%     SL = 2.5; 
%     QL =[rhoL*SL; rhoL*uL*SL; eL*SL] 
%     PL = ((gam-1)/SL)*(QL(3)-0.5*(QL(2)^2)/QL(1));
%     cL = sqrt(gam*PL/rhoL); 
%     uR = uRo/cLo; 
%     rhoR = rhoRo/rhoLo; 
%     eR = (rhoRo*ERo)/(rhoLo*cLo^2);
%     SR = 1.5; 
%     QR =[rhoR*SR; rhoR*uR*SR; eR*SR] 
%     PR = ((gam-1)/SR)*(QR(3)-0.5*(QR(2)^2)/QR(1));
%     cR = sqrt(gam*PR/rhoR); 
    
elseif (bcFlag == 2)
    
    uLo = 82.6934; 
    PLo = 96084.92;
    rhoLo = 1.1287;
    ELo = PLo/(gam-1) + 0.5*uLo^2;  
    cLo = sqrt(gam*PLo/rhoLo);
    SLo = 2.5; 
    QLo =[rhoLo*SLo; rhoLo*uLo*SLo; rhoLo*ELo*SLo]; 
    uRo = 151.6196; 
    PRo = 84973.94; 
    rhoRo = 1.0538; 
    ERo = PRo/(gam-1) + 0.5*uRo^2;  
    cRo = sqrt(gam*PRo/rhoRo); 
    SRo = 1.5; 
    QRo =[rhoRo*SRo; rhoRo*uRo*SRo; rhoRo*ERo*SRo];  

    % Non-dimensionalized 
    uL = uLo/cLo; 
    rhoL = rhoLo/rhoLo; 
    eL = (rhoLo*ELo)/(rhoLo*cLo^2);
    SL = 2.5; 
    QL =[rhoL*SL; rhoL*uL*SL; eL*SL]; 
    PL = ((gam-1)/SL)*(QL(3)-0.5*(QL(2)^2)/QL(1));
    cL = sqrt(gam*PL/rhoL); 
    uR = uRo/cLo; 
    rhoR = rhoRo/rhoLo;
    eR = (rhoRo*ERo)/(rhoLo*cLo^2);
    SR = 1.5; 
    QR =[rhoR*SR; rhoR*uR*SR; eR*SR]; 
    PR = ((gam-1)/SR)*(QR(3)-0.5*(QR(2)^2)/QR(1));
    cR = sqrt(gam*PR/rhoR); 
    
elseif (bcFlag == 3) 
    
    uLo=0;
    PLo = 10^5; 
    rhoLo = 1.0;
    ELo = PLo/(gam-1) + 0.5*uLo^2; 
    cLo = sqrt(gam*PLo/rhoLo);
    SLo = 2.5;
    QLo =[rhoLo*SLo; rhoLo*uLo*SLo; rhoLo*ELo*SLo]; 
    uRo=0; 
    PRo = 10^4; 
    rhoRo = 0.125;
    ERo = PRo/(gam-1) + 0.5*uRo^2; 
    cRo = sqrt(gam*PRo/rhoRo); 
    SRo = 2.5; 
    QRo =[rhoRo*SRo; rhoRo*uRo*SRo; rhoRo*ERo*SRo];  

    % Non-dimensionalized 
    uL = uLo/cLo;
    rhoL = rhoLo/rhoLo; 
    eL = (rhoLo*ELo)/(rhoLo*cLo^2);
    SL = 2.5;
    QL =[rhoL*SL; rhoL*uL*SL; eL*SL]; 
    PL = ((gam-1)/SL)*(QL(3)-0.5*(QL(2)^2)/QL(1));
    cL = sqrt(gam*PL/rhoL); 
    uR = uRo/cLo; 
    rhoR = rhoRo/rhoLo; 
    eR = (rhoRo*ERo)/(rhoLo*cLo^2);
    SR = 2.5; 
    QR =[rhoR*SR; rhoR*uR*SR; eR*SR];
    PR = ((gam-1)/SR)*(QR(3)-0.5*(QR(2)^2)/QR(1));
    cR = sqrt(gam*PR/rhoR); 
    
end

% Convert conservative flow variables back to dimensioned quantities
for j=1:nodes
    rho(j) = Q(3*j-2)*rhoLo/S(j);
    u(j) = (Q(3*j-1)/Q(3*j-2))*cLo;
    e(j) = (Q(3*j)/S(j))*rhoLo*cLo^2;
    
    Q(3*j-2) = rho(j)*S(j);
    Q(3*j-1) = rho(j)*u(j)*S(j);
    Q(3*j) = e(j)*S(j);
    
    P_d(j) = ((gam-1)/S(j))*(Q(3*j)-0.5*(Q(3*j-1)^2)/Q(3*j-2));
    rho_d(j) = Q(3*j-2)/S(j);
    u_d(j) = Q(3*j-1)/Q(3*j-2);
    c_d(j) = sqrt(gam*P_d(j)*S(j)/Q(3*j-2));
    M_d(j) = u_d(j)/c_d(j);
    T_d(j) = P_d(j)/(R*rho_d(j));   
    e_d(j) = Q(3*j)/S(j);    
end
