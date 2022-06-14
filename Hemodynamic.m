%% [Hemodynamic equation] 
% This script make a diferentian Hemodynamic 
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


function [p,q] = Hemodynamic(Z, U, P_SD, A,step)

        %In this implementation, the number of cortical regions (nRegions)
        %is inferred by size of Conectivity matrix (A)
        nRegions  = size(A,1);  
        % This is the number of timesteps of temporal input series (experimental inputs U ).
        simulationLength = size(U, 2);
        %TODO: Ok, use the raw dimention (2) works now but, it can be more
        %elegant.
        params = 6;   %Parametros 
        x = zeros(nRegions,params);
         
        %Priority values (#Regiones, numero params)
        kj = [0.69;0.63];
        yj = [0.28;0.34];
        tj = [2.11;1.97];
        tv = [4.12;0.92];
        
        %% Vasodilatory Signal  ----------------------
        s = zeros(nRegions,simulationLength-1);
        spriori = random('norm',0,0.000001,nRegions,1);
        s = [spriori,s];
        
        %% Inflow ------------------------------------
        f_in =  zeros(nRegions,simulationLength-1);
        f_in_priori = random('norm',0,0.000001,nRegions,1);
        f_in = [f_in_priori,f_in,];
        
        %% Outflow -----------------------------------
        f_out =  zeros(nRegions,simulationLength-1);
        f_out_priori = random('norm',0,0.000001,nRegions,1);
        f_out = [f_out_priori,f_out,];
 
        %% 
        tp =  zeros(nRegions,simulationLength-1);
        tp_priori = random('norm',0,0.000001,nRegions,1);
        tp = [tp_priori,tp,];
        
        tq =  zeros(nRegions,simulationLength-1);
        tq_priori = random('norm',0,0.000001,nRegions,1);
        tq = [tq_priori,tq];
        
        th =  zeros(nRegions,simulationLength-1);
        th_priori = random('norm',0,0.000001,nRegions,1);
        th = [th_priori,th];
        
        tvBallon =  zeros(nRegions,simulationLength-1);
        tvBallon_priori = random('norm',0,0.000001,nRegions,1);
        tvBallon = [tvBallon_priori,tvBallon];
        
        E =  zeros(nRegions,simulationLength-1);
        E_priori = random('norm',0,0.000001,nRegions,1);
        E = [E_priori,E];
        
        H = [0.64 0.32 2.00 0.32 0.32 2.00];
                    
        % signal decay
        sd = 0.64*exp(kj);
        % autoregulatory feedback
        af = 0.32.*exp(yj);
        % transit time
        tt = 2.00.*exp(tj);
        % viscoelastic time constant 
        tv = 2.00.*exp(tv); 
        
        for t = 2:simulationLength  
             % implement differential state equation f = dx/dt (hemodynamic)
             %--------------------------------------------------------------------------
             %Neuronal activity in source region j, zj causes an increase
             %in vasodilatoria signal sj.
             s_dot = Z(:,t-1) - sd.*s(:,t-1) - af.*(f_in(:,t-1)-0);             %Ec. 2.1 Tak  x(:,2)
             %El flujo de entrada responde en proporcion a esto ....
             f_in_dot = s(:,t-1);                                               %Ec. 2.2 Tak  x(:,3)
             
             %The rate of blood volume v change as... (TV)
             tvBallon_dot=(f_in(:,t-1)- (tvBallon(:,t-1)));                      %Ec. 3.2 Tak  x(:,4)
             
            
             E_f = 1-((0.5).^(f_in(:,t-1)));                                      %Ec. 4.2 Tak  no state
             
             %deoxyHb(tq) & totalHb(tp)
             tq_dot   = ((f_in(:,t-1).*(E(:,t-1)/0.32))-((tvBallon(:,t-1)).*(tq(:,t-1))./tvBallon(:,t-1)));
             tp_dot  = (f_in(:,t-1) - tvBallon(:,t-1)  .* tp(:,t-1)./tvBallon(:,t-1));         % x(:,6)


             %Euler
             s(:,t) = s(:,t-1) + step * s_dot;                         % x1
             f_in(:,t) = f_in(:,t-1) + step * f_in_dot;                % x2
             tvBallon(:,t) = tvBallon(:,t-1) + step * tvBallon_dot;    % x4
             E(:,t) = E_f;                                             % x5
             tq(:,t) = tq(:,t-1) + step * tq_dot;                      % x6
             tp(:,t) = tp(:,t-1) + step * tp_dot;                      % x3


        end
         
         q=f_in;
         plot(transpose(tq));
         shg;
         p=s;
         q=tq;
         
         
end
