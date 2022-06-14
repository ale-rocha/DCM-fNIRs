%% [Unit Test DCM-Taks] 
% This script make a replic of Taks modulatory example.
%% Framework
% 
% 
%% By:
% Instituto Nacional de Astrofísica Óptica y Electrónica
% Departamento de ciencias computacioanles.
% A.Rocha-Solache F.Orihuela-Espina, G.Rodríguez-Gómez
% rochasolache@inaoep.mx
%% Log:
% December, 2021
%   
%%
%% Display configurations
name_example = "Tak - unit test";   %Only a experiment name to display
verbose_plot = true;                %Plot results?
%% Model Tak params
%freq = params_series();
[A,B,C]= params_theta("Tak");
[U, timestamps] = getinputs(freq, 5, 25, 2); 

%% DCM start!
%% Neurodynamics
[Z] = Neurodynamics(A,B,C,U, 1/freq);

%% Hemodynamic
P_SD = [0.5 0.5 0.5 3];
[P,Q] = Hemo(Z, U, P_SD, A,1/freq);

%% Optic
Noise = 0;
[OR] = OpticLib( P,Q,U,A,Noise); 


%% Display results
if (verbose_plot == true)
    BilinearPlotThetaA(A,B,C,U,Z,P,OR);
end
  
