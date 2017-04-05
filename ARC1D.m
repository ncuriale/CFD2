% AER 1318
% Nathan Curiale
% 1002506781

clear all
clc
format long g

% Define constants
gam=1.4;% Specific heat ratio
R=287; % Specific gas constant for air
L=10; % Length of nozzle in meters
nodes=49; % Number of interior mesh points
dx=L/(nodes+1); % node spacing
tol=1e-14; % Residual convergence tolerance
pause_time=0.02; % Pause time for plot during convergence

% Enter flow case to be solved
bcFlag = input('Specify the flow problem to be solved (1-3): ');

% Enter method for solving system
diagFlag = input('Specify non-diagonal or diagonal form (1-2): ');

% Enter method for solving system
diffFlag = input('Specify Implicit Euler or 2nd Order Backward (1-2): ');

% Include boundary conditions
[uL, PL, rhoL, eL, SL, QL, cL, uR, PR, rhoR, eR, SR, QR, cR]= BCs(bcFlag);

% Calculate values of nozzle cross sectional area 
[S, S_p] = nozzleArea(nodes, dx, bcFlag);

% Calculate initial value of Q
Q = zeros(3*nodes,1);
Qiter = zeros(3*nodes,1);
for j=1:(nodes+1)/2
    Q(3*j-2) = rhoL*S(j);
    Q(3*j-1) = rhoL*uL*S(j);
    Q(3*j) = eL*S(j);
end
for j=(nodes+1)/2:nodes
    if (bcFlag<=2)
        Q(3*j-2) = rhoL*S(j);
        Q(3*j-1) = rhoL*uL*S(j);
        Q(3*j) = eL*S(j);
    elseif (bcFlag==3)
        Q(3*j-2) = rhoR*S(j);
        Q(3*j-1) = rhoR*uR*S(j);
        Q(3*j) = eL*S(j);
    end
end
Qiter(:,1) = Q; Qiter(:,2) = Q;%holds Qn and Qn-1 for 2nd order backwards

% Flow parameters at intial Q
[P, rho, u, M, T, c, e] = flowParam(S, Q, nodes);

%Set parameters based on case selected
if (bcFlag == 1)
    CFL = 100000;
    dt_nd = 0;   
elseif (bcFlag == 2)
    CFL = 40;
    dt_nd = 0;
elseif (bcFlag == 3)
    maxspeed=0;
    for j=1:nodes
        if ((abs(u(j))+c(j))>maxspeed)
            maxspeed=abs(u(j))+c(j);
        end
    end
    time = input('Time(s) for shock-tube solution: ');
    CFL = 1;
    dt_nd = CFL*dx/maxspeed;    % Time step in seconds
    time_steps = time/dt_nd;
    time_steps = round(time_steps); % Round number of time steps 
end

% Implicit euler to find solution
iter=0;
if (diagFlag == 1)
    if (bcFlag <= 2) %nozzle flow problems
        tic %loop timer
        for n=1:1000
            [LHS] = buildLHS(Q, S, S_p, nodes, dx, CFL, gam, bcFlag, dt_nd, diagFlag, diffFlag);
            [RHS] = buildRHS(Q, S, S_p, nodes, dx, CFL, gam, bcFlag, dt_nd, diagFlag, diffFlag, Qiter);
            normRHS(n) = norm(RHS);%store convergence history
            [iter norm(RHS) normRHS(n)/normRHS(1)]%show convergence
            
            if (normRHS(n)/normRHS(1)<tol)%check if converged to tolerance
                break;
            end
            
            dQ = LHS\RHS;%advance solution
            if (diffFlag == 2 && n==1)
                Qiter(:,1) = dQ + Q;%Use result of dQ as implicit euler and move into Qn 
            elseif (diffFlag == 2)
                Qiter(:,2) = Qiter(:,1);%Move Qn of this iteration to Qn-1
                Qiter(:,1) = dQ + Q;%Move Qn+1 of this iteration to Qn
            end
            Q = dQ + Q;%update Q
            iter = iter+1;
            
        end
        toc %end loop timer
    else %shock-tube problem
        tic %loop timer
        for n=1:time_steps
            [LHS] = buildLHS(Q, S, S_p, nodes, dx, CFL, gam, bcFlag, dt_nd, diagFlag, diffFlag);
            [RHS] = buildRHS(Q, S, S_p, nodes, dx, CFL, gam, bcFlag, dt_nd, diagFlag, diffFlag, Qiter);
            [iter norm(RHS)]%show convergence
            normRHS(n) = norm(RHS);%store convergence history
            
            dQ = LHS\RHS;%advance solution
            if (diffFlag == 2 && n==1)
                Qiter(:,1) = dQ + Q;%Use result of dQ as implicit euler and move into Qn 
            elseif (diffFlag == 2)
                Qiter(:,2) = Qiter(:,1);%Move Qn of this iteration to Qn-1
                Qiter(:,1) = dQ + Q;%Move Qn+1 of this iteration to Qn
            end
            Q = dQ + Q;%update Q
            iter = iter+1;

        end
        toc %end loop timer
     end
elseif (diagFlag == 2)
    if bcFlag <= 2 %nozzle flow problems
        tic %loop timer
        for n=1:1000     
            [LHS] = buildLHS(Q, S, S_p, nodes, dx, CFL, gam, bcFlag, dt_nd, diagFlag, diffFlag);
            [RHS] = buildRHS(Q, S, S_p, nodes, dx, CFL, gam, bcFlag, dt_nd, diagFlag, diffFlag, Qiter);
            normRHS(n) = norm(RHS);%store convergence history
            [iter norm(RHS) normRHS(n)/normRHS(1)]%show convergence
            
            if (normRHS(n)/normRHS(1)<tol)%check if converged to tolerance
                break;
            end
            
            Z = LHS\RHS;%advance solution
            
            % Calculate flow parameters at current time step
            [P, rho, u, M, T, c, e] = flowParam(S, Q, nodes);
            for j=1:nodes
                [X lambda Xinv] = diagFluxJac(j, P, rho, u, c, e, gam, S, S_p, nodes, bcFlag);
                dQ(3*j-2:3*j,1) = X*Z(3*j-2:3*j,1);
            end
            
            if (diffFlag == 2 && n==1)
                Qiter(:,1) = dQ + Q;%Use result of dQ as implicit euler and move into Qn 
            elseif (diffFlag == 2)
                Qiter(:,2) = Qiter(:,1);%Move Qn of this iteration to Qn-1
                Qiter(:,1) = dQ + Q;%Move Qn+1 of this iteration to Qn
            end
            Q = dQ + Q;%update Q
            iter = iter+1;
            
        end
        toc %end loop timer
    else %shock-tube problem
        tic %loop timer
        for n=1:time_steps
            [LHS] = buildLHS(Q, S, S_p, nodes, dx, CFL, gam, bcFlag, dt_nd, diagFlag, diffFlag);
            [RHS] = buildRHS(Q, S, S_p, nodes, dx, CFL, gam, bcFlag, dt_nd, diagFlag, diffFlag, Qiter);
            
            [iter norm(RHS)] %show convergence   
            Z = LHS\RHS;%advance soultion
            
            % Calculate flow parameters at current time step
            [P, rho, u, M, T, c, e] = flowParam(S, Q, nodes);
            
            for j=1:nodes
                [X lambda Xinv] = diagFluxJac(j, P, rho, u, c, e, gam, S, S_p, nodes, bcFlag);
                dQ(3*j-2:3*j,1) = X*Z(3*j-2:3*j,1);
            end
            
            if (diffFlag == 2 && n==1)
                Qiter(:,1) = dQ + Q;%Use result of dQ as implicit euler and move into Qn 
            elseif (diffFlag == 2)
                Qiter(:,2) = Qiter(:,1);%Move Qn of this iteration to Qn-1
                Qiter(:,1) = dQ + Q;%Move Qn+1 of this iteration to Qn
            end
            Q = dQ + Q;%update Q
            iter = iter+1;
                        
        end
        toc %end loop timer
    end    
end 

fprintf('Total Iterations: %4f\n', iter);

% Final flow solution quantities
[P, rho, u, M, T, c, e] = flowParam(S, Q, nodes);

% Calculate exact solutions
if (bcFlag==1)   
    [U_exact, RHO_exact, P_exact, M_exact] = nozzleExact(L, dx);    
elseif (bcFlag == 2)    
    [U_exact, RHO_exact, P_exact, M_exact] = transnozzleExact(L, dx);       
elseif (bcFlag == 3)   
    [U_exact, RHO_exact, P_exact, M_exact] = shocktubeExact(L, time, dx);      
end

% Plot numerical and exact solutions
plotSol(P, M, rho, P_exact, M_exact, RHO_exact, bcFlag, nodes)

% Plot convergence
if (bcFlag<=2)
    plotConv(normRHS, bcFlag, iter)
end

% Error between exact and numerical solutions
M_error = zeros(nodes,1);
if (bcFlag<=2)
    P_error = zeros(nodes,1);
else
    Rho_error = zeros(nodes,1);
end
for j=1:nodes
    M_error(j) = M_error(j) + abs(M(j)-M_exact(j+1));
    if (bcFlag<=2)
        P_error(j) = P_error(j) + abs(P(j)-P_exact(j+1));
    else
        Rho_error(j) = Rho_error(j) + abs(rho(j)-RHO_exact(j+1));
    end
end
M_errorNorm = norm(M_error);
M_errorSum = sum(M_error);
if (bcFlag<=2)
    P_errorNorm = norm(P_error);
    P_errorSum = sum(P_error);
    fprintf('Mach Error Norm: %4f\n', M_errorNorm);
    fprintf('Mach Error Sum: %4f\n', M_errorSum);
    fprintf('Pressure Error Norm: %4f\n', P_errorNorm);
    fprintf('Pressure Error Sum: %4f\n', P_errorSum);
else
    Rho_errorNorm = norm(Rho_error);
    Rho_errorSum = sum(Rho_error);
    fprintf('Time Step Size: %4f\n', dt_nd);
    fprintf('Mach Error Norm: %4f\n', M_errorNorm);
    fprintf('Mach Error Sum: %4f\n', M_errorSum);
    fprintf('Density Error Norm: %4f\n', Rho_errorNorm );
    fprintf('Density Error Sum: %4f\n', Rho_errorSum );
end





