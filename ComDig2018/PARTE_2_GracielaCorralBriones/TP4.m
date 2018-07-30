%--------------------------------------------------------------------------
%   Simulador de un Sistema de Comunicación con Modulación 2PAM en un canal|
%                       con ruido gaussiano. Version-1.0                   |
%              Laboratorio de Comunicaciones Digitales FCEFYN-UNC          |
%                                                                          |
%                                    By Martin (martin.ayarde@unc.edu.ar)  |
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Transmisor                                              |Canal
% |-------|   |---------|   |---------|   |-----------|  |  |-----------|
% |Gen. de|   |Conver.  |   | Sobre-  |   |  Filtro   |  |  | Canal con |
% |de bits|-->|de Bits a|-->|muestreo |-->|Conformador|--|->|  ruido    |-->
% |       |   |símbolos |   |         |   |           |  |  | Gaussiano |
% |-------|   |---------|   |---------|   |-----------|  |  |-----------|
%-------------------------------------------------------------------------
%Receptor                       |------------------------> Análisis SER Vs EsNo
%    |-----------|   |--------| | |-----------|
%    |  Filtro   |   |  Sub-  | | |Conver.    |
% -->| Apareado  |-->|muestreo|-->|de símbolos|----------> Análisis BER Vs EbNo
%    |           |   |        |   |a bits     |
%    |-----------|   |--------|   |-----------|
%--------------------------------------------------------------------------
clear all;close all;clc;
%Parámetros
Channel	=   'H'; % rcos, H, 613, A, B
% 'OFF, 'ZeroISI' ''MMSE' 'Adapt' 
% Para canal 6.13 DFE y Adapt solo
EQU='MMSE'; 
K = 6; % K de los eq
DataLength = 10000; %4000000;         % Cantidad de símbolos a transmitir
FreqSymbol = 1.2e3;          % Frecuencia de transmisión de símbolos
% Cantidad de ceros entre símbolos para asegurar cero ISI(interferencia entre símbolos)
% Modificar su valor entre 6 y 12 y justificar degradación respecto a la curva teórica.
%Evaluar NumZeros en 19 y justificar la mejora en la curva de SER obtenida
if strcmp(Channel,'rcos')
    NumZeros = 9;                
elseif strcmp(Channel,'H')
    %Parametros de H
    NumZeros = 1; % Muestrado al doble de la tasa de simbolo                
else
    %Parametros de H
    NumZeros = 0; % Muestrado al doble de la tasa de simbolo                
end
SwitchChannel = 'ON';        % ON: Enciende el canal de ruido gaussiano,
% OFF: Apaga el canal de ruido gaussiano
EsNodB = 2:10;               % Determina la relación en decibeles entre la energia del símbolo Es y la densidad espectral de potencia de ruido No
EsNo = 10.^(EsNodB/10);
SwitchSignal = 'RandomData';  % 'RandomData': Datos binarios generados con una distribición uniforme
% 'ProbeData': Datos binarios de ceros y unos intercalados
% 'RealData': Datos binarios generados a partir de una imagen
RxFilter = 'Matched';         %'Matched':filtro receptor que se encuentra apareado con el filtro transmisor
%'NotMatched':filtro receptor que no esta apareado con el filtro transmisor
ImageFile = 'Image.jpg';      %Cargar imagen a transmitir
%--------------------------------------------------------------------------
% Diseño del filtro conformador de la señal
RollOff = 0.50;
Span = 10;  % Esto es el span del raised cos en duracion de symbolo
Sps = 10;   % Muestras por span, o sea, cada 10 (Sps) debe ir un simnbolo
% O sea que el largo de el pulso es Span*Sps=100, 50 para cada lado del
% central
Shape = 'sqrt'; % Raiz del coseno realsado

SRRCFilter = rcosdesign(RollOff,Span,Sps,Shape);  % Diseña un filtro del tipo raiz coseno realzado
RCFilter = rcosdesign(RollOff,Span,Sps,'normal'); % Diseña un filtro del tipo coseno realzado
%SRRCFilter = rcosine(1,Sps,Shape,RollOff,floor(Span/2));  para versiones de
%matlab que no soporten rcosdesign reemplazar la linea anterior por esta linea

if length(EsNodB)== 1
    %figure(1)
    %subplot 211
    %stem([0:length(SRRCFilter)-1]/(FreqSymbol*Sps),SRRCFilter);grid;
    %[ssf,FFT] = plotspect(1/(FreqSymbol*Sps),SRRCFilter);
    %subplot 212
    %plot(ssf,20*log10(FFT))
    %grid on;
    %xlabel('Frequency [Hz]'); ylabel('Magnitude')
    % Graficos de cascada de filtros tx rx
    figure(1)
    subplot 211
    stem([0:length(RCFilter)-1]/(FreqSymbol*Sps),RCFilter);grid;
    [ssf,FFT] = plotspect(1/(FreqSymbol*Sps),RCFilter);
    subplot 212
    plot(ssf,20*log10(FFT))
    grid on;
    xlabel('Frequency [Hz]'); ylabel('Magnitude')
end

%--------------------------------------------------------------------------
if strcmp(SwitchSignal,'RealData')
    Image = Image2Bits(ImageFile); %Convierte una imagen jpg en un vector de bits
    %pause(5)
end
NumError = zeros(length(DataLength),1);
for k=1:length(EsNodB)
    %--------------------------------------------------------------------------
    % Transmisor digital:
    %rng(1024)   %Utilizado para generar siempre la misma secuencia de bits aleatorios
    if strcmp(SwitchSignal,'ProbeData')
        SymbolsAux = ones(1,DataLength);
        SymbolsAux(2:2:end) = -1;
    elseif strcmp(SwitchSignal,'RandomData')
        SymbolsAux = 2*(rand(1,DataLength)>0.5)-1;  % Generación de símbolos +/-1
    elseif strcmp(SwitchSignal,'RealData')
        SymbolsAux = 2*(Image>0.5)-1;
    else
        display('Opcion no valida')
    end
    if length(EsNodB)== 1
        figure(2)
        stem(SymbolsAux,'r'), grid;
    end
    Symbols = upsample(SymbolsAux,NumZeros+1);  % Agregado de ceros entre los símbolos
    % Filtro cascada rcos
    
    if strcmp(Channel,'rcos')
        TxPAMSignal = conv(Symbols,SRRCFilter);       % Conformación de la señal a transmitir
        Maximo = max(TxPAMSignal);
        if length(EsNodB)== 1
            figure(3)
            subplot 211
            plot(TxPAMSignal);grid;
            subplot 212
            [pxx w] = pwelch(TxPAMSignal,1024,512,1024*2);
            plot(w*FreqSymbol*Sps/(2*pi),10*log10(sqrt(pxx)));grid on
            xlabel('Frecuencia [Hz]')
            ylabel('PSD [dB]');
            title(' Densidad Espectral de Potencia de |p[n]|^2')

        end
        %--------------------------------------------------------------------------
        %Canal Gausiano:
        if strcmp(SwitchChannel,'ON')
            Noise = randn(1,length(TxPAMSignal))./sqrt(2*EsNo(k));
        elseif strcmp(SwitchChannel,'OFF')
            Noise = zeros(1,length(TxPAMSignal));
        else
            display('Opción no valida')
        end
        NoisySignal = TxPAMSignal+Noise;
        %--------------------------------------------------------------------------
        %Receptor Digital:
        if strcmp(RxFilter,'Matched')
            RxPAMSignal = conv(NoisySignal,SRRCFilter);
        elseif strcmp(RxFilter,'NotMatched')
            RxPAMSignal = conv(NoisySignal,RCFilter);
        end
        if length(EsNodB)==1
            figure(4)
            stem(RxPAMSignal);grid;
        end
    elseif strcmp(Channel,'H')
        % ----------------------
        % Canal H con distorcion
        % ----------------------
        
        %   Ajustes del canal
        T  = NumZeros + 1;  %Simbolos cada 2
        Ts = T/2; 
        L = 5; % Extremos
        t = -L*T:Ts:L*T;
        H = 1./(1+((T/2)*t).^2);
        % Respuesta al impulso del canal
         figure
         stem(H)
        
        % Senal a la salida del Filtro Receptor
        TxPAMSignal = conv(Symbols,H);
        % Ruido de la senal
        if strcmp(SwitchChannel,'ON')
            Noise = randn(1,length(TxPAMSignal))./sqrt(2*EsNo(k));
        elseif strcmp(SwitchChannel,'OFF')
            Noise = zeros(1,length(TxPAMSignal));
        else
            display('Opción no valida') 
        end
        RxPAMSignal = TxPAMSignal+Noise;
        % constelacion de simbolos recibidos.   
        %ri =  11 ; % Retardo hasta primer simbolo L*T+1 % 1 para adp
        %cSymb = 1000;
        %figure('Name',['Constelacion para ' int2str(EsNo(k))])
        %stem(RxPAMSignal(ri:T:cSymb),zeros(1,length(ri:T:cSymb)),'b:*')
       
     
   elseif strcmp(Channel,'613')
        % ----------------------
        % Canal 613 con distorcion
        % ----------------------
        
        %   Ajustes del canal
        X = [0.05 -0.063 0.088 -0.126 -0.25 0.9047 0.25 0 0.126 0.038 0.088];
        T=1;
        %figure;
        %stem(X);
        
        % Senal a la salida del Filtro Receptor
        TxPAMSignal = conv(SymbolsAux,X);
        % Ruido de la senal
        if strcmp(SwitchChannel,'ON')
            Noise = randn(1,length(TxPAMSignal))./sqrt(2*EsNo(k));
        elseif strcmp(SwitchChannel,'OFF')
            Noise = zeros(1,length(TxPAMSignal));
        else
            display('Opción no valida') 
        end
        RxPAMSignal = TxPAMSignal+Noise;
        % constelacion de simbolos recibidos.   
        %ri =  11 ; % Retardo hasta primer simbolo L*T+1 % 1 para adp
        %cSymb = 1000;
        %figure('Name',['Constelacion para ' int2str(EsNodB(k))])
        %stem(RxPAMSignal(ri:T:cSymb),zeros(1,length(ri:T:cSymb)),'b:*')
       
     elseif strcmp(Channel,'A')
        % ----------------------
        % Canal A con distorcion
        % ----------------------
        
        %   Ajustes del canal
        X = [0.04 -0.05 0.07 -0.21 -0.5 0.72 0.63 0 0.21 0.03 0.07];
        T=1;
        %figure('Name','Respuesta del canal A');
        %stem(X);
        
        % Senal a la salida del Filtro Receptor
        TxPAMSignal = conv(SymbolsAux,X);
        % Ruido de la senal
        if strcmp(SwitchChannel,'ON')
            Noise = randn(1,length(TxPAMSignal))./sqrt(2*EsNo(k));
        elseif strcmp(SwitchChannel,'OFF')
            Noise = zeros(1,length(TxPAMSignal));
        else
            display('Opción no valida') 
        end
        RxPAMSignal = TxPAMSignal+Noise;
        % constelacion de simbolos recibidos.   
        %ri =  1 ;
        %cSymb = 200;
        %figure('Name',['Constelacion para ' int2str(EsNodB(k))])
        %stem(RxPAMSignal(ri:T:cSymb),zeros(1,length(ri:T:cSymb)),'b:*')
    elseif strcmp(Channel,'B')
        % ----------------------
        % Canal B con distorcion
        % ----------------------
        
        %   Ajustes del canal
        X = [0.407 0.815 0.407];
        T=1;
        %figure('Name','Respuesta del canal B');
        %stem(X);
        
        % Senal a la salida del Filtro Receptor
        TxPAMSignal = conv(SymbolsAux,X);
        % Ruido de la senal
        if strcmp(SwitchChannel,'ON')
            Noise = randn(1,length(TxPAMSignal))./sqrt(2*EsNo(k));
        elseif strcmp(SwitchChannel,'OFF')
            Noise = zeros(1,length(TxPAMSignal));
        else
            display('Opción no valida') 
        end
        RxPAMSignal = TxPAMSignal+Noise;
        % constelacion de simbolos recibidos.   
        %ri =  1 ;
        %cSymb = 500;
        %figure('Name',['Constelacion para ' int2str(EsNodB(k))])
        %stem(RxPAMSignal(ri:T:cSymb),zeros(1,length(ri:T:cSymb)),'b:*')
       
     
    end
    % -----------------------------------------
    %   Estapa de equalizacio para el canal H
    % -----------------------------------------
    if strcmp(Channel,'H') 
        % Armo matrix X
        MX = zeros((2*K+1),(2*K+1));
        for m = -K:1:K
            for n = -K:1:K
                MX(m+K+1,n+K+1) = (1/(1+(2*(m*T-n*Ts)/T)^2));
            end
        end
        
        if strcmp(EQU,'ZeroISI') 
            qSol  = zeros(2*K+1,1);
            qSol(K+1,1) = 1;
            c_opt = inv(MX)*qSol
            % Imprimimos respuesta del EQ y respuesta a la salida del EQ
             figure
             stem(c_opt);
             figure
            outEQ = conv(c_opt,H);
            % Downsampling
            stem(outEQ(2:2:end));
            title('ecualizado');
            outFilter = conv(c_opt,RxPAMSignal);
            RxPAMSignal = outFilter(K+1:end);
        end
        
        if strcmp(EQU,'MMSE') 
            No = 1/EsNo(k);
            Ray  = MX(K+1,:);
            Ry   = MX'*MX + (No/2)*eye(2*K+1);
            c_opt = inv(Ry)*Ray';
            % Imprimimos respuesta del EQ y respuesta a la salida del EQ
            %figure
            %stem(c_opt);
            % Downsampling
            figure('Name',['Pulso para  ' int2str(EsNo(k))])            %outEQ = conv(c_opt,H);
            outEQ = conv(c_opt,H);
            stem(outEQ(1:2:end));
            outFilter = conv(c_opt,RxPAMSignal);
            RxPAMSignal = outFilter(K+1:end);
        end
        if strcmp(EQU,'Adapt') 
            delta  = 1/10000;
            
            Nent   = length(Symbols);
            
            retardo = (length(H)-(2*K+1))/2;
            if retardo>=0
                RxPAMSignal = RxPAMSignal(retardo + 1 :end);
            else
                RxPAMSignal = [zeros(1,-retardo) RxPAMSignal];
            end
            
            % Etapa de entrenamiento con los mismos datos de la img
            z_k = zeros(1,Nent);
            c_est = zeros(1,2*K+1);
            c_est(1,K+1) = 1;
            c_est_history = zeros(Nent-2*K,2*K+1);
            
            for i=1:2:Nent-2*K
                y_k = RxPAMSignal(i:i+2*K);
                z_k(i) = c_est*y_k.';
                e_k = Symbols(i) - z_k(i);
                c_est_history(i,:) = [c_est];
                c_est = c_est + delta*e_k*y_k;
            end
            
            figure('Name',['Evolucion de coef Para ' int2str(EsNodB(k))])
            for i=1:2*K+1
                subplot(K+1,2,i) 
                plot(c_est_history([1:10:100000],i))
                coef = i-(K+1);
                title(['Coeficiente ' int2str(-coef)])

            end
            % Fin etapa de entrenamiento
            
            for i=1:2:Nent-2*K
                y_k = RxPAMSignal(i:i+2*K);
                z_k(i) = c_est*y_k.';
                a_k = (z_k(i)>0)*2-1;
                e_k = a_k - z_k(i);
                c_est = c_est + delta*e_k*y_k;
            end
            
            RxPAMSignal = z_k;

            
            % Imprimimos respuesta del EQ y respuesta a la salida del EQ
        	figure('Name',['Para ' int2str(EsNodB(k))])
            subplot(1,2,1) 
            stem(flip(c_est));
            title('Coeficientes')
            outEQ = conv(H, [flip(c_est) zeros(1,retardo)]);
            subplot(1,2,2) 
            stem(outEQ(1:1:end));
            title('Pulso equalizado')
            
        end
        
        
        % constelacion de simbolos recibidos.   
        ri =  11 ; % Retardo hasta primer simbolo L*T+1 % 1 para adp
        cSymb = 1000;
        figure('Name',['Constelacion para ' int2str(EsNodB(k))])
        stem(RxPAMSignal(ri:T:cSymb),zeros(1,length(ri:T:cSymb)),'b:*')
        
    end

    
    % -----------------------------------------
    %   Estapa de equalizacio para el canal 613 A y B
    % -----------------------------------------
    if strcmp(Channel,'613') % 613 A B OFF
        if strcmp(EQU,'Adapt') 
            delta  = 1/1000;
            Nent   = length(Symbols);
           
            retardo = (length(X)-(2*K+1))/2;
            if retardo>=0
                RxPAMSignal = RxPAMSignal(retardo + 1 :end);
            else
                RxPAMSignal = [zeros(1,-retardo) RxPAMSignal];
            end
            
            % Etapa de entrenamiento con los mismos datos de la img
            z_k = zeros(1,Nent);
            c_est = zeros(1,2*K+1);
            c_est(1,K+1) = 1;
            c_est_history = zeros(Nent-2*K,2*K+1);
            
            for i=1:1:Nent-2*K
                y_k = RxPAMSignal(i:i+2*K);
                z_k(i) = c_est*y_k.';
                e_k = Symbols(i) - z_k(i);
                c_est_history(i,:) = [c_est];
                c_est = c_est + delta*e_k*y_k;
            end
        	figure('Name',['Evolucion de coef Para ' int2str(EsNodB(k))])
            for i=1:2*K+1
                subplot(K+1,2,i) 
                plot(c_est_history([1:10:end],i))
                coef = i-(K+1);
                title(['Coeficiente ' int2str(-coef)]); % - por el flip

            end

            
            % Fin etapa de entrenamiento
            
            for i=1:1:Nent-2*K
                y_k = RxPAMSignal(i:i+2*K);
                z_k(i) = c_est*y_k.';
                a_k = (z_k(i)>0)*2-1;
                e_k = a_k - z_k(i);
                c_est = c_est + delta*e_k*y_k;
            end
            
            RxPAMSignal = z_k;
            
            % Imprimimos respuesta del EQ y respuesta a la salida del EQ
        	figure('Name',['Para ' int2str(EsNodB(k))])
            subplot(1,2,1) 
            stem(flip(c_est));        	
            title('Coeficientes')
            outEQ = conv(X, [flip(c_est) zeros(1, retardo)]); % Ajustar los ceros segun K
            subplot(1,2,2) 
            stem(outEQ(1:1:end));
            title('Pulso equalizado')


        elseif strcmp(EQU,'DFE') 
            display(['llegue a ' int2str(EsNodB(k))])
            delta  = 1/10000;
            %delta  = 0.07169;
            Nent   = length(Symbols);
            K1 = 6;
            K2 = 1;

           retardo = (length(X)-(2*K1+1))/2;
           if retardo>=0
                RxPAMSignal = RxPAMSignal(retardo + 1 :end);
            else
                RxPAMSignal = [zeros(1,-retardo) RxPAMSignal];
            end
            
            % Etapa de entrenamiento con los mismos datos
            
            % FF filter
            c_est_ff = zeros(1,2*K1+1);
            c_est_ff(1,K1+1) = 1;
            c_est_history_ff = zeros(Nent-2*K1,2*K1+1);
            % FB filter
            c_est_fb = zeros(1,2*K2+1);
            c_est_fb(1,2*K2+1) = 1;
            c_est_history_fb = zeros(Nent-2*K1,2*K2+1);
            % secuencia detectasda
            a_k = zeros(1,Nent);
            z_k = zeros(1,Nent);
            
            for i=1:1:Nent-2*K1
                % Entrada al FF filter y salida
                y_k_ff = RxPAMSignal(i:i+2*K1);
                z_k_ff = c_est_ff*y_k_ff.';
                % Entrada el FB filter y salida
                if i==1
                    a_k_fb = zeros(1,2*K2+1);
                elseif i-1<(2*K2+1) % i-1 cantidad de simbolos detectados
                    a_k_fb = [zeros(1,2*K2+1-(i-1)) a_k(1:i-1)];
                else
                    a_k_fb = a_k(i-(2*K2+1):i-1);
                end
                z_k_fb = c_est_fb*a_k_fb.';
                % Retroalimentacion
                z_k(i) = z_k_ff - z_k_fb;
                % Detector
                a_k(i) = (z_k(i)>0)*2-1;
                e_k = Symbols(i) - z_k(i);
                %Calculo de coef
                c_est_ff = c_est_ff + delta*e_k*y_k_ff;
                c_est_fb = c_est_fb - delta*e_k*a_k_fb*10;
                c_est_history_ff(i,:) = [c_est_ff];
                c_est_history_fb(i,:) = [c_est_fb];
            end
            
        	figure('Name',['Evolucion de coef Para FF' int2str(EsNodB(k))])
            for j=1:2*K1+1
                subplot(K1+1,2,j) 
                plot(c_est_history_ff([1:100:end],j))
                coef = j-(K1+1);
                title(['Coeficiente ' int2str(-coef)])
            end
            figure('Name',['Evolucion de coef Para FB' int2str(EsNodB(k))])
            for j=1:2*K2+1
                subplot(K2+1,2,j) 
                plot(c_est_history_fb([1:100:end],j))
                coef = j-(K2+1);
               title(['Coeficiente ' int2str(-coef)])
            end
            z_k = zeros(1,Nent);

            % Fin etapa de entrenamiento
            
            for i=1:1:Nent-2*K1
                % Entrada al FF filter y salida
                y_k_ff = RxPAMSignal(i:i+2*K1);
                z_k_ff = c_est_ff*y_k_ff.';
                % Entrada el FB filter y salida
                if i==1
                    a_k_fb = zeros(1,2*K2+1);
                elseif i-1<(2*K2+1) % i-1 cantidad de simbolos detectados
                    a_k_fb = [zeros(1,2*K2+1-(i-1)) a_k(1:i-1)];
                else
                    a_k_fb = a_k(i-(2*K2+1):i-1);
                end
                z_k_fb = c_est_fb*a_k_fb.';
                % Retroalimentacion
                z_k(i) = z_k_ff - z_k_fb;
                % Detector
                a_k(i) = (z_k(i)>0)*2-1;
                e_k = a_k(i) - z_k(i);
                %Calculo de coef
                c_est_ff = c_est_ff + delta*e_k*y_k_ff;
                c_est_fb = c_est_fb - delta*e_k*a_k_fb*10;
                c_est_history_ff(i,:) = [c_est_ff];
                c_est_history_fb(i,:) = [c_est_fb];
            end
            
            
            RxPAMSignal = z_k;
            
            % Imprimimos respuesta del EQ y respuesta a la salida del EQ
        	%figure('Name',['Para ' int2str(EsNodB(k))])
            %subplot(1,2,1) 
            %stem(flip(c_est));        	
            %title('Coeficientes')
            %outEQ = conv(X, [flip(c_est) zeros(1, retardo)]);
            % Downsampling
            %subplot(1,2,2) 
            %plot(outEQ(1:1:end));
            %title('Pulso equalizado')
            %figure
            %subplot(1,2,2) 
            %outEQ = conv(X, [c_est_ff]);
            %plot(outEQ(1:1:end));
            %figure
            

        end
        
        % constelacion de simbolos recibidos.   
        ri =  1 ; % Retardo hasta primer simbolo L*T+1 % 1 para adp
        cSymb = 200;
        figure('Name',['Constelacion para ' int2str(EsNodB(k))])
        stem(RxPAMSignal(ri:T:cSymb),zeros(1,length(ri:T:cSymb)),'b:*')
        
    end
    
    
    
    
    
    
    
    
    
    
    if strcmp(Channel,'rcos')
        % Retardo para tener en cuenta los transitorios producidos por las 
        % conv al inicio y al final de la transmisión
        RetardoIni = 101; % floor((Sps*Span*2)/2+1) Por 2 de c/filtro /2+1 centrarlo en mx
    else
        if strcmp(EQU,'Adapt') 
            RetardoIni = 1;
        elseif strcmp(EQU,'DFE') 
            RetardoIni = 1;
        elseif strcmp(EQU,'OFF') 
            RetardoIni = 1;
        else
            RetardoIni = 11; % El canal son 21 muestras, centrando el pulso de max E 10+1
        end
    end

          
    ErrorRetardoIni = 0;   % Modificar su valor entre -9 y 9 y justificar degradación respecto a la curva teórica
    ErrorFrecSimbol = 0;   % Modificar su valor entre 0 y -3 y justificar degradación respecto a la curva teórica
    Detector = 2*(RxPAMSignal(RetardoIni+ErrorRetardoIni:NumZeros+1+ErrorFrecSimbol:end)>0)-1; % convierte las muestras de los símbolos recibidas en valores +-1
    Detector = Detector(1:length(SymbolsAux));
    if length(EsNodB)==1
        figure(5)
        stem(Detector);grid;
    end
    NumError(k) = sum(SymbolsAux~=Detector); %Contador de errores
    Rxbits = Detector>0;  %Conversor de símbolos a bits
    if strcmp(SwitchSignal,'RealData')
        fHat = Bits2Image(ImageFile,Rxbits,num2str(EsNodB(k)));
    end
end
pause(2)
%--------------------------------------------------------------------------
%Estimación del tiempo utilizado para la transmisión de la información
M = 2;     %Tipo de modulacion 2PAM, 4PAM, ect
TotalTime = length(Rxbits)/(FreqSymbol*log2(M)*60);
%--------------------------------------------------------------------------
%Análisis del desempeño del sistema de comunicaciones mediante curvas de Bit Error Rate (BER):
simSER = NumError/length(SymbolsAux);  %Estimación del BER obtenido
figure
semilogy(EsNodB,0.5*erfc(sqrt(EsNo)),'k');
hold on; semilogy(EsNodB,simSER,'b:o')
title('Curva de tasa de error de Símbolos (SER)')
xlabel('Es/No (dB)')
ylabel('Probabilidad de error')
grid on

