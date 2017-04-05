function [AOut, BOut, COut, DOut] = calcOutBlocks(j, dx, dt, Q, nodes, S, S_p, gam, bcFlag)   

% Flow parameters at current Q
[P, rho, u, M, T, c, e] = flowParam(S, Q, nodes);

% Dissipation terms at node j
[c4, c2] = dissTerms(j, u, c, P, dx, nodes, S, gam, bcFlag);

% Calculate the A outflow 
AOut = (c4(1))*eye(3,3);
  
% Calculate the B outflow
[fluxJac, sourceJac] = calcJacs((j-1), P, rho, u, c, e, gam, S, S_p, nodes, bcFlag);
BOut = (-1/(2*dx))*fluxJac - (c4(2) + 3*c4(1) + c2(1))*eye(3,3);
       
% Calculate the C outflow
[fluxJac, sourceJac] = calcJacs(j, P, rho, u, c, e, gam, S, S_p, nodes, bcFlag);
if (j==nodes)
    COut = (1/dt)*eye(3,3) - sourceJac + (2*c4(2) + 3*c4(1) + c2(2) + c2(1))*eye(3,3);
else
    COut = (1/dt)*eye(3,3) - sourceJac + (3*c4(2) + 3*c4(1) + c2(2) + c2(1))*eye(3,3);
end

% Calculate the D outflow
[fluxJac, sourceJac] = calcJacs((j+1), P, rho, u, c, e, gam, S, S_p, nodes, bcFlag);
if (j==nodes)
    DOut=0; 
else
    DOut = (1/(2*dx))*fluxJac - (3*c4(2) + c4(1) + c2(2))*eye(3,3);   
end    



        
    




