function [RHS] = buildRHS(Q, S, S_p, nodes, dx, CFL, gam, bcFlag, dt_nd, diagFlag, diffFlag, Qiter)
% RHS of the implicit Euler time marching method for quasi-1D Euler equations

% Flow parameters at current Q
[P, rho, u, M, T, c, e] = flowParam(S, Q, nodes);

% time step
dt=zeros(nodes,1);
for j=1:nodes 
    if bcFlag <= 2;
        dt(j) = CFL*dx/(abs(u(j))+c(j));       
    else 
        dt(j) = dt_nd;
    end
end
dt=min(dt);

RHS=zeros(3*nodes,1);
for j=1:nodes
    [G] = calcSourceTerm(P, S_p, j);%source term at node j
    [FD] = calcFluxDiff(j, Q, S, nodes, gam, dx, bcFlag);%2nd order differenced flux at node j
    [D] = calcDissipation(j, Q, S, gam, nodes, dx, bcFlag);%dissipation at node j
    
    %2nd Order Backwards Difference
    if (diffFlag==2)
        A=(1/(2*dt))*(Qiter((3*j-2:3*j),1)-Qiter((3*j-2:3*j),2));
    else
        A=zeros(3,1);
    end
  
    if (diagFlag==1)
        RHS_jsum = A + G - FD + D;%calculate the RHS of Implicit Euler
    elseif (diagFlag==2)
        %inverse of flux jacobian eigenvalue matrix
        [X lambda Xinv] = diagFluxJac(j, P, rho, u, c, e, gam, S, S_p, nodes, bcFlag);
        RHS_jsum = A + G - FD + D;%calculate the RHS of Implicit Euler
        RHS_jsum = Xinv*RHS_jsum;
    end
    
    RHS((3*j-2:3*j),1) = RHS_jsum;% Store result of RHS 
    
end

   