function [X lambda Xinv] = diagFluxJac(k, P, rho, u, c, e, gam, S, S_p, nodes, bcFlag);

% Include boundary conditions
[uL, PL, rhoL, eL, SL, QL, cL, uR, PR, rhoR, eR, SR, QR, cR] = BCs(bcFlag);
         
if (k == 0)

    X(1,1) = 1; 
    X(2,1) = (uL+cL); 
    X(3,1) = (uL^2/2 + cL*uL + cL^2/(gam-1)); 
    X(1,2) = 1; 
    X(2,2) = (uL-cL); 
    X(3,2) = (uL^2/2 - cL*uL + cL^2/(gam-1));
    X(1,3) = 1;
    X(2,3) = uL;
    X(3,3) = uL^2/2;
    
    lambda(1,1) = uL+cL; 
    lambda(1,2) = 0; 
    lambda(1,3) = 0;
    lambda(2,1) = 0; 
    lambda(2,2) = uL-cL;
    lambda(2,3) = 0;
    lambda(3,1) = 0; 
    lambda(3,2) = 0;
    lambda(3,3) = uL;
    
    Xinv=inv(X);
    
elseif (k <= nodes)

    X(1,1) = 1; 
    X(2,1) = (u(k)+c(k)); 
    X(3,1) = (u(k)^2/2 + c(k)*u(k) + c(k)^2/(gam-1)); 
    X(1,2) = 1;
    X(2,2) = (u(k)-c(k)); 
    X(3,2) = (u(k)^2/2 - c(k)*u(k) + c(k)^2/(gam-1));
    X(1,3) = 1;
    X(2,3) = u(k);
    X(3,3) = u(k)^2/2;
    
    lambda(1,1) = u(k)+c(k); 
    lambda(1,2) = 0; 
    lambda(1,3) = 0;
    lambda(2,1) = 0; 
    lambda(2,2) = u(k)-c(k); 
    lambda(2,3) = 0;
    lambda(3,1) = 0; 
    lambda(3,2) = 0;
    lambda(3,3) = u(k);
    
    Xinv=inv(X); 
    
else 
    
    X(1,1) = 1; 
    X(2,1) = (uR+cR); 
    X(3,1) = (uR^2/2 + cR*uR + cR^2/(gam-1)); 
    X(1,2) = 1; 
    X(2,2) = (uR-cR); 
    X(3,2) = (uR^2/2 - cR*uR + cR^2/(gam-1));
    X(1,3) = 1;
    X(2,3) = uR;
    X(3,3) = uR^2/2;
    
    lambda(1,1) = uR+cR; 
    lambda(1,2) = 0; 
    lambda(1,3) = 0;
    lambda(2,1) = 0; 
    lambda(2,2) = uR-cR; 
    lambda(2,3) = 0;
    lambda(3,1) = 0; 
    lambda(3,2) = 0;
    lambda(3,3) = uR;
    
    Xinv=inv(X); 
    
end
          



