function [fluxJac, sourceJac] = calcJacs(k, P, rho, u, c, e, gam, S, S_p, nodes, bcFlag);

% Include boundary conditions
[uL, PL, rhoL, eL, SL, QL, cL, uR, PR, rhoR, eR, SR, QR, cR] = BCs(bcFlag);

% Calculate the flux Jacobian matrix
if (k==0)    
    fluxJac(1,1) = 0;
    fluxJac(1,2) = 1;
    fluxJac(1,3) = 0;
    fluxJac(2,1) = ((gam-3)/2)*uL^2;
    fluxJac(2,2) = -(gam-3)*uL;
    fluxJac(2,3) = (gam-1);
    fluxJac(3,1) = -gam*eL*uL/rhoL + (gam-1)*uL^3;
    fluxJac(3,2) = gam*eL/rhoL - 1.5*(gam-1)*uL^2;
    fluxJac(3,3) = gam*uL;    
elseif (k<=nodes)
    fluxJac(1,1) = 0;
    fluxJac(1,2) = 1;
    fluxJac(1,3) = 0;
    fluxJac(2,1) = ((gam-3)/2)*u(k)^2;
    fluxJac(2,2) = -(gam-3)*u(k);
    fluxJac(2,3) = (gam-1);
    fluxJac(3,1) = -gam*e(k)*u(k)/rho(k) + (gam-1)*u(k)^3;
    fluxJac(3,2) = gam*e(k)/rho(k) - 1.5*(gam-1)*u(k)^2;
    fluxJac(3,3) = gam*u(k);    
else
    fluxJac(1,1) = 0;
    fluxJac(1,2) = 1;
    fluxJac(1,3) = 0;
    fluxJac(2,1) = ((gam-3)/2)*uR^2;
    fluxJac(2,2) = -(gam-3)*uR;
    fluxJac(2,3) = (gam-1);
    fluxJac(3,1) = -gam*eR*uR/rhoR + (gam-1)*uR^3;
    fluxJac(3,2) = gam*eR/rhoR - 1.5*(gam-1)*uR^2;
    fluxJac(3,3) = gam*uR;       
end

% Calculate the source term Jacobian matrix
if (k==0)  
    sourceJac(1,1) = 0;
    sourceJac(1,2) = 0;
    sourceJac(1,3) = 0;
    sourceJac(2,1) = ((gam-1)/2)*uL^2;
    sourceJac(2,2) = -(gam-1)*uL;
    sourceJac(2,3) = (gam-1);
    sourceJac(3,1) = 0;
    sourceJac(3,2) = 0;
    sourceJac(3,3) = 0;    
elseif (k<=nodes)
    sourceJac(1,1) = 0;
    sourceJac(1,2) = 0;
    sourceJac(1,3) = 0;
    sourceJac(2,1) = ((gam-1)/2)*u(k)^2;
    sourceJac(2,2) = -(gam-1)*u(k);
    sourceJac(2,3) = (gam-1);
    sourceJac(3,1) = 0;
    sourceJac(3,2) = 0;
    sourceJac(3,3) = 0;    
else
    sourceJac(1,1) = 0;
    sourceJac(1,2) = 0;
    sourceJac(1,3) = 0;
    sourceJac(2,1) = ((gam-1)/2)*uR^2;
    sourceJac(2,2) = -(gam-1)*uR;
    sourceJac(2,3) = (gam-1);
    sourceJac(3,1) = 0;
    sourceJac(3,2) = 0;
    sourceJac(3,3) = 0;     
end

% If solving for the steady-state nozzle flow cases
if bcFlag <= 2
    if (k==0)
        sourceJac = (1/SL)*(-0.6)*sourceJac;
    elseif (k<=nodes)
        sourceJac = (1/S(k))*(S_p(k))*sourceJac;
    else
        sourceJac = (1/SR)*(0.2)*sourceJac;   
    end    
else    
    if (k==0)
        sourceJac = (1/SL)*(0)*sourceJac;
    elseif (k<=nodes)
        sourceJac = (1/S(k))*(S_p(k))*sourceJac;
    else
        sourceJac = (1/SR)*(0)*sourceJac;   
    end      
end





