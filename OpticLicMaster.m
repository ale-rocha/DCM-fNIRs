%% [Optic equation] 
% This script make a diferentian optic 
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

function [response_temporal,respuesta] = OpticLicMaster(P,Q,A)
         

         simulationLength = size(P, 2); % Longjtud de la simulacion
         nRegions = size(P, 1);
         fuentes_opticas = 6;
         %In principle, all optical channels are sensitive to all neuronal
         %sources. In practice, the sensitivity of channels to sources is
         %governed by the sensitivity matrix. It provides -N channels-
         %measurements of optical density changesat time t and wavelegth
         %(lamda1) to (lambda2).
         wavelengths = 2;
         
         disp(" [Optic] Dimensiones de la simulacion ");
         disp(size(P));
         disp(" [Optic] Longitudes de onda ");
         disp(wavelengths);

         
         %% optics parameters
         %--------------------------------------------------------------------------
         %   N(1) - oxygen saturation                                            SO2
         %   N(2) - baseline total-hb concentration [mM]                          P0
         %   N(3) - cortical weighting factor                                      w
         %--------------------------------------------------------------------------
         N = [0.65 0.071 2]; 
         % Normalized to their baseline concentration
         base_hbr = N(2) .* (1-N(1));
         
         %% calculate changes in hemoglobin concentration  Ec. 7  Tak------
         DQ = base_hbr.*(Q-1); % dxy-Hb
         DP = N(2).*(P-1); % total Hb
         DH = DP - DQ; % oxy-Hb
         disp(" [Optic] Dimensiones incremntos de deoxy ");
         disp(size(DH));
        
         % parameter for correction of pial vein contamination ---- W ----
         W = N(3).* exp([0.5.^.76;0.7.^.1]);        % (2X2)
         W = kron(ones(1,2), W); 
         epsilon=[4.26,4.26;3.88,3.88];
         epsilon=[4.26;3.38];
         disp("Dimensiones de la matriz de contaminacion de pial venas");
         disp(size(W));
         
         %TODO : fuentes distribuidas
  
         
         %S = ones(fuentes_opticas,nRegions);
         S = random('norm',5,3,fuentes_opticas,nRegions); %Aleatori samples
         
         disp("Dimensiones de la matriz de sensibilidad");
         disp(size(S));
         
        % calculate optical density changes  para cada tiempo t y para cada
        % longitud de onda
        response_temporal = ones(wavelengths,fuentes_opticas,simulationLength);
        for t = 1:simulationLength
            disp("Computando el tiempo-----------------------------------");
            disp(t);
            response_lambda = ones(wavelengths,fuentes_opticas,1);
            for i = 1:wavelengths
                %Temporalmente asignare unos valores aleatorios de la matri
                %de sensibilidad.
                
                %g(:,i) = (pv .* (S * [sh(i,:) sq(i,:)])) * M(i,:);
                incrementos = [DH(i,t),DH(i,t); DQ(i,t),DQ(i,t)];
                disp("Dimension de los incrementos");
                disp(size(incrementos));
                temp = S * incrementos;   % ------- 1
                %      (2x2)*(2X1) = (2X1)
                disp("Dimension de tenp");
                disp(size(temp));
                disp("Dimension de W");
                disp(size(W));
                temp = temp*W ;
                %      (2x2)*(2X1) = (2X1)
                disp("Dimension epsilon");
                disp(size(epsilon));
                disp("Temp");
                disp(size(temp));
                temp=(temp*epsilon);

                response_lambda(i,:,:) = temp; % respuesta lambda 
            end
            response_temporal (1,:,t)= response_lambda(1,:,:); %Guardando longitus de onda 1
            response_temporal (2,:,t)= response_lambda(2,:,:); %Guardando longitus de onda 2
        end
        
        respuesta.lambda1=response_temporal(1,:,:);
        respuesta.lambda2=response_temporal(2,:,:);
        
        disp("Respuesta temporal");
        disp(size(response_temporal));
        response_temporal = reshape(response_temporal,[wavelengths*fuentes_opticas,simulationLength]);
        response_temporal = transpose(response_temporal);
        plot(response_temporal);
        shg;
        
        %esta función regresa un array de [receptores X tiempo]
        %cada receptor observa en 2 longitudes de onda.
        %
        %  reciber1 -lambda1 -------------------------> time
        %  reciber1 -lambda2
        %  reciber2 -lambda1
        %  reciber2 -lambda2
        %   ....
end
