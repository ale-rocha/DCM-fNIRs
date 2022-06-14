%% [Neurodynamic equation] 
% This script make a diferentian neurodynamic calculation
%% Framework
% 
% 
%% Autors:
% Instituto Nacional de Astrofísica Óptica y Electrónica
% Departamento de ciencias computacioanles.
% A.Rocha-Solache F.Orihuela-Espina, G.Rodríguez-Gómez
% rochasolache@inaoep.mx
%% Log activity:
% 14 - Aug - 2021 : Creation file
% 1 - Feb - 2022 : Resolving to do task
%   
%% Biblio
% [Tak S.] - Tak,S., Kempny,A., Friston,K.J., Leff,A.P., & Penny,W.D. 
%            (2015). Dynamic causal modelling for functional near-infrared
%            spectroscopy. Neuroimage, 111, 338-349.

function [Z] = Neurodynamics(A, B, C, U, step)
   
   % TODO: Implement the legacy file to check concordance among matrix
   %size and cortical regions(nRegions) and experimental inputs(M) desirable.
   
   %In this implementation, the number of cortical regions (nRegions)
   %is inferred by size of Conectivity matrix (A)
   nRegions  = size(A,1);
   %Same case to M (experimental inputs), the size is inferred from
   %modulatory matrix U size.
   M = size(U,1); 
   
   % This is the number of timesteps of temporal input series (experimental inputs U ).
   simulationLength = size(U,2);
   %TODO: Ok, use the raw dimention (2) works now but, it can be more
   % elegant.
 
   % Neurual activity matriz
   Z = zeros(nRegions,simulationLength-1);
   % We set de priori values in 0's to all cortical areas in neural
   % activity matriz (Z). A gaussian priori stimation can result in a nicely gap.
   Zpriori = zeros(nRegions,1);
   Z = [Z,Zpriori]; 

   %This is a very simple implementation of euler method
   %I know, can be more precise
   for timestep = 2:simulationLength  %% For each timestep...
       
       % I use the variable J to express the effect of experimental input
       % along the modulatory matriz in all cortical regions.
       J = zeros(nRegions); %Expected shape (nRegions x nRegions)
       for corticalReg1 = 1:M %% along cortical regions regions (M) row
           for corticalReg2 = 1:M %% along cortical regions regions (M) columms
                if (corticalReg2 == corticalReg2)
                    % "... To enforce constraints on the parameters being estimated, we make
                    % use of latent variable ... J_(self_conexions), are  strictly negative with a typical
                    % value of 0.5" [Tak .S page:340]
                    J(corticalReg1,corticalReg2) = -0.5 * exp ( (U(corticalReg1,timestep) * B(corticalReg1,corticalReg2)));
                else
                    %"... For the off-diagonal terms, we have no such constraints.."
                    % [Tak .S page:340]
                    J(corticalReg1,corticalReg2) =  exp ( (U(corticalReg1,timestep) * B(corticalReg1,corticalReg2)));
                end
           end
       end
       
       %Self-inhibition is SI
       SI = diag(A); 
       A =  A - diag(exp(SI)/2 + SI);
      
       %Finally compute Zdot (neural activity)!
       Zdot = (A + J) * Z(:,timestep-1) + C*U(:,timestep-1); 
       
       %Euler iteration
	   Z(:,timestep) = Z(:,timestep-1) + step * Zdot;
       
       %Finaly, Dynamic B changes for each U acording to a gaussian
       %distribution
       % Mean and std form Tak S. pag. 
       priorMeanB = 0; priorStdB = 0.2;
       B = normrnd(priorMeanB,priorStdB,[nRegions,nRegions]);
   end   
  
   
end
