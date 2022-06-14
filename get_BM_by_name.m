%% [get_schema_bm] 
% This script return the outputs for a bilinear model schema

%% Autors:
% Instituto Nacional de Astrofísica Óptica y Electrónica
% Departamento de ciencias computacioanles.
% A.Rocha-Solache F.Orihuela-Espina, G.Rodríguez-Gómez
% rochasolache@inaoep.mx
%% Log activity:
% 31 - May - 2021 : Creation file
%   
%% Biblio
% [Tak S.] - Tak,S., Kempny,A., Friston,K.J., Leff,A.P., & Penny,W.D. 
%            (2015). Dynamic causal modelling for functional near-infrared
%            spectroscopy. Neuroimage, 111, 338-349.


function [SMA,M1,status,log] = get_BM_by_name(name,verbose)

    [params_series,params_dcm,status_name,v] = legacy_name_cdm_model(name,verbose);
    
    if status_name 
        
        %% Model Tak params
        freq = params_series.Freq;     
        
        %% Generate Series
        [U, timestamps] = getinputs(freq, 5, 25, 2); 

        %% DCM start!
        %% Neurodynamics
        [Z] = Neurodynamics(params_dcm.A,params_dcm.B,params_dcm.C,U, 1/freq);

        %% Hemodynamic
        P_SD = [0.5 0.5 0.5 3];
        [P,Q] = Hemodynamic(Z, U, P_SD, params_dcm.A,1/freq);

        %% Optic
        [SMA,M1] = OpticLib( P,Q,U,params_dcm.A,params_series.Noise); 

        %% Display results
        %if (verbose_plot == true)
        %    BilinearPlotThetaA(params_dcm.A,params_dcm.B,params_dcm.C,U,Z,P,OR);
        %end
        status = true;
        log = "Generado con exito";
    else
        log = "Incorrect name provided";
        status = false;
    end
end
