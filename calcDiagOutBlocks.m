function [AOut, BOut, COut, DOut] = calcDiagOutBlocks(j, dx, dt, Q, nodes, S, S_p, gam, bcFlag)                                       

% Flow parameters at current Q
[P, rho, u, M, T, c, e] = flowParam(S, Q, nodes);

% Dissipation terms at node j
[c4, c2] = dissTerms(j, u, c, P, dx, nodes, S, gam, bcFlag);

% Calculate A outflow
AOut = (c4(1))*eye(3,3);
  
% Calculate B outflow
[X lambda Xinv] = diagFluxJac((j-1), P, rho, u, c, e, gam, S, S_p, nodes, bcFlag);
BOut = (-1/(2*dx))*lambda - ((c4(2) + 3*c4(1) + c2(1)))*eye(3,3);
       
% Calculate C outflow 
if (j==nodes)
    COut = (1/dt)*eye(3,3) + (2*c4(2) + 3*c4(1) + c2(2) + c2(1))*eye(3,3);
else
    COut = (1/dt)*eye(3,3) + (3*c4(2) + 3*c4(1) + c2(2) + c2(1))*eye(3,3);
end

% Calculate the D outflow
[X lambda Xinv] = diagFluxJac((j+1), P, rho, u, c, e, gam, S, S_p, nodes, bcFlag);
if (j==nodes)
    DOut = 0;
else
    DOut = (1/(2*dx))*lambda - ((3*c4(2) + c4(1) + c2(2)))*eye(3,3);
end    



        
    




