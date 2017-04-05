function [S, S_p] = nozzleArea(nodes, dx, bcFlag)
%   Defines cross-sectional area along tube depending on case selected.
%
%   Steady-state nozzle flow cases use a converging/diverging
%   Unsteady shock tube flow cases use a tube of constant cross-sectional area

S = zeros(nodes,1); % Holds nozzle cross sectional values

if (bcFlag <=2)
    for j = 1:(5/dx)
        S(j) = 1 + 1.5*(1-(j*dx)/5)^2;
    end
    for j = ((5/dx)+1):nodes
        S(j) = 1 + 0.5*(1-(j*dx)/5)^2;
    end

    S_p = zeros(nodes,1); % Holds nozzle cross sectional derivatives (dS/dx)

    for j = 1:(5/dx)
        S_p(j) = -3/5 + 3*j*dx/25;
    end
    for j = ((5/dx)+1):nodes
        S_p(j) = -1/5 + j*dx/25;
    end
    
else 
    S = 2.5*ones(nodes,1);   
    S_p = zeros(nodes,1);       
end


 