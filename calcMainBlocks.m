function [A, B, C, D, E] = calcMainBlocks(j, dx, dt, Q, nodes, S, S_p, gam, bcFlag)

% Flow parameters at current Q
[P, rho, u, M, T, c, e] = flowParam(S, Q, nodes);

% Dissipation terms at node j
[c4, c2] = dissTerms(j, u, c, P, dx, nodes, S, gam, bcFlag);

% Calculate A
A = (c4(1))*eye(3,3);
    
% Calculate B
[fluxJac, sourceJac] = calcJacs((j-1), P, rho, u, c, e, gam, S, S_p, nodes, bcFlag);
B = (-1/(2*dx))*fluxJac - ((c4(2) + 3*c4(1) + c2(1)))*eye(3,3);
    
% Calculate C
[fluxJac, sourceJac] = calcJacs(j, P, rho, u, c, e, gam, S, S_p, nodes, bcFlag);
C = (1/dt)*eye(3,3) - sourceJac + ((3*c4(2) + 3*c4(1) + c2(2) + c2(1)))*eye(3,3);       

% Calculate D
[fluxJac, sourceJac] = calcJacs((j+1), P, rho, u, c, e, gam, S, S_p, nodes, bcFlag);
D = (1/(2*dx))*fluxJac - ((3*c4(2) + c4(1) + c2(2)))*eye(3,3);

% Calculate E
E = (c4(2))*eye(3,3);

        
    




