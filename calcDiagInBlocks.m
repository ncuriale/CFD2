function [BIn, CIn, DIn, EIn] = calcDiagInBlocks(j, dx, dt, Q, nodes, S, S_p, gam, bcFlag)

% Flow parameters at current  Q
[P, rho, u, M, T, c, e] = flowParam(S, Q, nodes);

% Dissipation terms at node j
[c4, c2] = dissTerms(j, u, c, P, dx, nodes, S, gam, bcFlag);

% Calculate B inflow
[X lambda Xinv] = diagFluxJac((j-1), P, rho, u, c, e, gam, S, S_p, nodes, bcFlag);                  
if (j==1)
    BIn = 0;
else
    BIn = (-1/(2*dx))*lambda - (c4(2) + 3*c4(1) + c2(1))*eye(3,3);
end
    
% Calculate C inflow
if (j==1)
    CIn = (1/dt)*eye(3,3) + ((3*c4(2) + 2*c4(1) + c2(2) + c2(1)))*eye(3,3);
else
    CIn = (1/dt)*eye(3,3) + ((3*c4(2) + 3*c4(1) + c2(2) + c2(1)))*eye(3,3);
end

% Calculate D inflow
[X lambda Xinv] = diagFluxJac((j+1), P, rho, u, c, e, gam, S, S_p, nodes, bcFlag);
DIn = (1/(2*dx))*lambda - (3*c4(2) + c4(1) + c2(2))*eye(3,3);

% Calculate E inflow
EIn = (c4(2))*eye(3,3);
        
    




