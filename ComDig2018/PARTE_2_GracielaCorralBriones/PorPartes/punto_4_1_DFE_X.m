% Ecualización adaptiva. Estudiar el desempeño de los siguientes 
% esquemas. Explicar los resultados a través de diagrama de 
% constelaciones, respuesta impulsiva y BER. Analizar la influencia
% del tamaño de los filtros adaptivos en el desempeño.
% Implementar los ecualizadores adaptivos lineal y realimentado 
% por decisiones para el canal del Ejemplo 6.12. Nota: el canal 
% discreto está muestreado a la velocidad del símbolo.
clear all;close all;clc;
%Parámetros
DataLength = 1000000;        % Cantidad de símbolos a transmitir
FreqSymbol = 1.2e3;         % Frecuencia de transmisión de símbolos
NumZeros = 0;               % Cantidad de ceros para sobremuestreo
SwitchChannel = 'ON';      % ON: Enciende el canal de ruido gaussiano,
                            % OFF: Apaga el canal de ruido gaussiano
EsNodB = 2:10;               % SNR en dB
EsNo = 10.^(EsNodB/10);
SwitchSignal = 'RandomData';% 'RandomData'
                            % 'ProbeData'
                            % 'RealData'
% K = 5;                      % 11 Coeficientes
K = 6;                      % 13 Coeficientes
%--------------------------------------------------------------------------
% Configuracion de graficos
Graficos = 'OFF';           % Apaga/enciende graficos de datos transmitidos
Constelation = 'OFF';
GraficosEqu = 'ON';          %Canal
SeleccionGrafico = 9;       % Determina grafico segun k
BERDiagram = 'ON';
EQU = 'ON';                 %Activa/Desactiva Ecualizador adaptivo
%--------------------------------------------------------------------------
NumError = zeros(length(DataLength),1);
for k=1:length(EsNodB)
    % Seleccion informacion a transmitir
    if strcmp(SwitchSignal,'ProbeData')
        SymbolsAux = ones(1,DataLength);
        SymbolsAux(2:2:end) = -1;
    elseif strcmp(SwitchSignal,'RandomData')
        SymbolsAux = 2*(rand(1,DataLength)>0.5)-1;
    elseif strcmp(SwitchSignal,'RealData')
        SymbolsAux = 2*(Image>0.5)-1;
    else
        display('Opcion no valida')
    end
    
    % Dibuja simbolos a transmitir sin sobremuestreo
    if(strcmp(Graficos,'ON')&& k == SeleccionGrafico)
        figure
        stem(SymbolsAux(1:20),'r'), grid;
        title('SIMBOLOS A TRANSMITIR, SIN SOBREMUESTREO');
    end
    
    % Sobremuestreo de simbolos
    Symbols = upsample(SymbolsAux,NumZeros+1);
    
    % Dibuja simbolos a transmitir con sobremuestreo
    if(strcmp(Graficos,'ON')&&k==SeleccionGrafico)
        figure
        stem(Symbols(1:200),'b');
        title('SIMBOLOS A TRANSMITIR, CON SOBREMUESTREO');
    end
    %----------------------------------------------------------------------
    % Canal X (613) con distorsion - DFE 
    %   Ajustes del canal
    X = [0.05 -0.063 0.088 -0.126 -0.25 0.9047 0.25 0 0.126 0.038 0.088];
    T=1;
    if(strcmp(GraficosEqu,'ON')&&k==SeleccionGrafico)
        figure
        stem(X);
        title('RESPUESTA AL IMPULSO DEL CANAL X');
    end

    % Senal a la salida del Filtro Receptor
    TxPAMSignal = conv(SymbolsAux,X);
    %----------------------------------------------------------------------
    %Canal Gausiano:
    if strcmp(SwitchChannel,'ON')
        Noise = randn(1,length(TxPAMSignal))./sqrt(2*EsNo(k));
    elseif strcmp(SwitchChannel,'OFF')
        Noise = zeros(1,length(TxPAMSignal));
    else
        display('Opción no valida')
    end
    %----------------------------------------------------------------------
    %Receptor Digital:
    RxPAMSignal = TxPAMSignal+Noise;
    
    if(strcmp(Graficos,'ON')&&k==SeleccionGrafico)
        figure
        stem(RxPAMSignal(1:200));grid;
        title('SENAL RECIBIDA DESPUES DEL FILTRO, CON SOBREMUESTREO');
    end  
    %----------------------------------------------------------------------
    % Ecualizador del canal X - DFE
     if strcmp(EQU,'ON') 
        display(['llegue a ' int2str(EsNodB(k))])
        delta  = 1/1000;
        %delta  = 0.07169;
        Nent   = length(Symbols);
        K1 = 6; %FF - Directo
        K2 = 1; %FB - Realimentado

%         retardo = (length(X)-(2*K1+1))/2;
%         
%            if retardo>=0
%                 RxPAMSignal = RxPAMSignal(retardo + 1 :end);
%             else
%                 RxPAMSignal = [zeros(1,-retardo) RxPAMSignal];
%            end
            
% Etapa de entrenamiento---------------------------------------------------

        % FF filter - Directo
        c_est_ff = zeros(1,2*K1+1);
        c_est_ff(1,K1+1) = 1;
        c_est_history_ff = zeros(Nent-2*K1,2*K1+1);
        % FB filter - Realimentado
        c_est_fb = zeros(1,2*K2+1);
        c_est_fb(1,K2+1) = 1;
        c_est_history_fb = zeros(Nent-2*K1,2*K2+1);
        % secuencia detectada
        a_k = zeros(1,Nent);
        z_k = zeros(1,Nent);
            
        for i=1:1:Nent-2*K1
            % Entrada al FF filter y salida
            y_k_ff = RxPAMSignal(i:i+2*K1); %Tomo tramos de la senal recibida
            z_k_ff = c_est_ff*y_k_ff.';     %C_estimado_FF * tramos de la senal recibida
            
            % Entrada el FB filter y salida
            if i==1
                a_k_fb = zeros(1,2*K2+1);
            elseif i-1<(2*K2+1) % i-1 cantidad de simbolos detectados
                a_k_fb = [zeros(1,2*K2+1-(i-1)) a_k(1:i-1)];
            else
                a_k_fb = a_k(i-(2*K2+1):i-1);
            end
            
            z_k_fb = c_est_fb*a_k_fb.';     %C_estimado_FB * tramos de la senal recibida
            % Retroalimentacion
            z_k(i) = z_k_ff - z_k_fb;
            % Detector
            a_k(i) = (z_k(i)>0)*2-1;    %Detector
            e_k = Symbols(i) - z_k(i); %Error 
            %Calculo de coef
            c_est_ff = c_est_ff + delta*e_k*y_k_ff; %Actualizacion coef Directo
            c_est_fb = c_est_fb - delta*e_k*a_k_fb*10; %Actualizacion coef Realimentado
            
            c_est_history_ff(i,:) = [c_est_ff];
            c_est_history_fb(i,:) = [c_est_fb];
        end
        
%Graficos de evolucion de los coeficientes---------------------------------
        if (mod(k,2) ~= 0) %Evolucion coeficientes para FF (Directo)  
            figure('Name',['Evolucion de coef Para FF' int2str(EsNodB(k))])
            for j=1:2*K1+1
                subplot(K1+1,2,j) 
                plot(c_est_history_ff([1:100:end],j))
                coef = j-(K1+1);
                title(['Coeficiente ' int2str(coef)])
            end
        end     
        
        if (mod(k,2) ~= 0) %Evolucion coeficientes para FB (Realimentado) 
            figure('Name',['Evolucion de coef Para FB' int2str(EsNodB(k))])
            for j=1:2*K2+1
                subplot(K2+1,2,j) 
                plot(c_est_history_fb([1:100:end],j))
                coef = j-(K2+1);
               title(['Coeficiente ' int2str(coef)])
            end
        end
%--------------------------------------------------------------------------
        z_k = zeros(1,Nent);

% Fin etapa de entrenamiento-----------------------------------------------
            
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
            e_k = a_k(i) - z_k(i); %Error (Salida de los filtros)
            %Calculo de coef
            c_est_ff = c_est_ff + delta*e_k*y_k_ff;
            c_est_fb = c_est_fb - delta*e_k*a_k_fb*10;
            
            c_est_history_ff(i,:) = [c_est_ff];
            c_est_history_fb(i,:) = [c_est_fb];
        end
%--------------------------------------------------------------------------            
            
        RxPAMSignal = z_k;

        % Imprimimos respuesta del EQ y respuesta a la salida del EQ
%         figure('Name',['Para ' int2str(EsNodB(k))])
%         subplot(1,2,1) 
%         stem(c_est);        	
%         title('Coeficientes')
%         outEQ = conv(X, c_est_fb);
%         % Downsampling
%         subplot(1,2,2) 
%         stem(outEQ(1:1:end));
%         figure
%         title('Pulso equalizado')
%         subplot(1,2,2) 
%         outEQ = conv(X, c_est_fb);
%         stem(outEQ(1:1:end));
        %figure
     end
     
     % constelacion de simbolos recibidos.   
        ri =  1 ; % Retardo hasta primer simbolo L*T+1 % 1 para adp
        cSymb = 200;
        figure('Name',['Constelacion para ' int2str(EsNodB(k))])
        stem(RxPAMSignal(ri:T:cSymb),zeros(1,length(ri:T:cSymb)),'b:*')
    
    %----------------------------------------------------------------------
    % Retardo para tener en cuenta los transitorios producidos por las 
    % conv al inicio y al final de la transmisión
    RetardoIni = 1; % floor((Sps*Span*2)/2+1) Por 2 de c/filtro /2+1 centrarlo en mx
    ErrorRetardoIni = 0;   % Modificar su valor entre -9 y 9 y justificar degradación respecto a la curva teórica
    ErrorFrecSimbol = 0;   % Modificar su valor entre 0 y -3 y justificar degradación respecto a la curva teórica
    Detector = 2*(RxPAMSignal(RetardoIni+ErrorRetardoIni:NumZeros+1+ErrorFrecSimbol:end)>0)-1; % convierte las muestras de los símbolos recibidas en valores +-1
    Detector = Detector(1:length(SymbolsAux));
    if(strcmp(Graficos,'ON')&&k==SeleccionGrafico)
        figure
        stem(Detector(1:20));grid;
        title('SENAL DETECTADA');
    end
    NumError(k) = sum(SymbolsAux~=Detector); %Contador de errores
    Rxbits = Detector>0;  %Conversor de símbolos a bits
    
    % Constelacion de simbolos recibidos.  
%     if (strcmp(Constelation,'ON')&&mod(k,2)~=0)
%         ri =  11 ; % Retardo hasta primer simbolo L*T+1 % 1 para adp
%         cSymb = 1000;
%         figure('Name',['Constelacion para ' int2str(EsNo(k))])
%         stem(RxPAMSignal(ri:T:cSymb),zeros(1,length(ri:T:cSymb)),'b:*')
%     end
end

%--------------------------------------------------------------------------
%Estimación del tiempo utilizado para la transmisión de la información
M = 2;     %Tipo de modulacion 2PAM, 4PAM, ect
TotalTime = length(Rxbits)/(FreqSymbol*log2(M)*60);
%--------------------------------------------------------------------------
%Análisis del desempeño del sistema de comunicaciones mediante curvas de Bit Error Rate (BER):
if strcmp(BERDiagram,'ON')
    simSER = NumError/length(SymbolsAux);  %Estimación del BER obtenido
    figure
    semilogy(EsNodB,0.5*erfc(sqrt(EsNo)),'k');
    hold on; semilogy(EsNodB,simSER,'b:o')
    title('Curva de tasa de error de Símbolos (SER)')
    xlabel('Es/No (dB)')
    ylabel('Probabilidad de error')
    grid on
end