function [LHS] = buildLHS(Q, S, S_p, nodes, dx, CFL, gam, bcFlag, dt_nd, diagFlag, diffFlag)
% LHS of the implicit Euler time marching method for quasi-1D Euler equations

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

%multiplier for 2nd order backwards scheme
if (diffFlag == 2)
    dt=2*dt/3;
end

for j=1:nodes
    % Calculate block matrices of LHS   
    if (j<=2)       
        if (diagFlag == 1)
            [BIn, CIn, DIn, EIn] = calcInBlocks(j, dx, dt, Q, nodes, S, S_p, gam, bcFlag);           
        elseif (diagFlag ==2)            
            [BIn, CIn, DIn, EIn] = calcDiagInBlocks(j, dx, dt, Q, nodes, S, S_p, gam, bcFlag);  
        end    
   elseif (j<(nodes-1))
        if (diagFlag == 1)
            [A, B, C, D, E] = calcMainBlocks(j, dx, dt, Q, nodes, S, S_p, gam, bcFlag);            
        elseif (diagFlag == 2)
            [A, B, C, D, E] = calcDiagMainBlocks(j, dx, dt, Q, nodes, S, S_p, gam, bcFlag);   
        end                  
    else       
        if (diagFlag == 1) 
            [AOut, BOut, COut, DOut] = calcOutBlocks(j, dx, dt, Q, nodes, S, S_p, gam, bcFlag);         
        elseif (diagFlag == 2)
            [AOut, BOut, COut, DOut] = calcDiagOutBlocks(j, dx, dt, Q, nodes, S, S_p, gam, bcFlag);           
        end                                    
    end

    
    % Insert block matrices into inflow  
    if (j==1)     
        LHS(1:3,1:3) = CIn;
        LHS(1:3,4:6) = DIn;
        LHS(1:3,7:9) = EIn;       
    elseif (j==2)       
        LHS(4:6,1:3) = BIn;
        LHS(4:6,4:6) = CIn;
        LHS(4:6,7:9) = DIn; 
        LHS(4:6,10:12) = EIn;
        
    % Insert block matrices into interior points      
    elseif (j<(nodes-1))        
        LHS(3*j-2:3*j,3*j-8:3*j-6) = A;
        LHS(3*j-2:3*j,3*j-5:3*j-3) = B;
        LHS(3*j-2:3*j,3*j-2:3*j) = C;
        LHS(3*j-2:3*j,3*j+1:3*j+3) = D; 
        LHS(3*j-2:3*j,3*j+4:3*j+6) = E; 
        
    % Insert block matrices into outflow    
    elseif (j==nodes-1)        
        LHS(3*j-2:3*j,3*j-8:3*j-6) = AOut;
        LHS(3*j-2:3*j,3*j-5:3*j-3) = BOut;
        LHS(3*j-2:3*j,3*j-2:3*j) = COut;
        LHS(3*j-2:3*j,3*j+1:3*j+3) = DOut;        
    elseif (j==nodes)        
        LHS(3*j-2:3*j,3*j-8:3*j-6) = AOut;
        LHS(3*j-2:3*j,3*j-5:3*j-3) = BOut;
        LHS(3*j-2:3*j,3*j-2:3*j) = COut;
        
    end
    
end



