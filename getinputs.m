%% Autors:
% Instituto Nacional de Astrofísica Óptica y Electrónica
% Departamento de ciencias computacioanles.
% A.Rocha-Solache F.Orihuela-Espina, G.Rodríguez-Gómez
% rochasolache@inaoep.mx
%% Log activity:
% 31-May-22 : Creation file
%
%   
%% Biblio
% [Tak S.] - Tak,S., Kempny,A., Friston,K.J., Leff,A.P., & Penny,W.D. 
%            (2015). Dynamic causal modelling for functional near-infrared
%            spectroscopy. Neuroimage, 111, 338-349.


function [U, timestamps] = BilinearModel_StimulusTrainGenerator(freq, action_time, rest_time, cycles)
    rest  = rest_time * freq;
    activation = action_time * freq;
    activation_ = activation;
    Time_period  = (action_time + rest_time) * cycles;
    simulationLength = Time_period * freq;
    U = zeros(2,simulationLength);
    
    lvel = 1;
    for i = 1:cycles
        for x = lvel:activation_
            U(:,x) = 1; 
            lvel = lvel + 1;
        end
        lvel = (activation+rest)*i;
        activation_ =  ((activation+rest)*i)+activation;
    end    
    
    timestamps = [0:1/freq:Time_period];
    timestamps(end)=[];
end