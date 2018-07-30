% Ecualizaci�n adaptiva. Estudiar el desempe�o de los siguientes 
% esquemas. Explicar los resultados a trav�s de diagrama de 
% constelaciones, respuesta impulsiva y BER. Analizar la influencia
% del tama�o de los filtros adaptivos en el desempe�o.
% Implementar los ecualizadores adaptivos lineal y realimentado 
% por decisiones para el canal del Ejemplo 6.12. Nota: el canal 
% discreto est� muestreado a la velocidad del s�mbolo.
clear all;close all;clc;
%Par�metros
DataLength = 100000;        % Cantidad de s�mbolos a transmitir
FreqSymbol = 1.2e3;       % Frecuencia de transmisi�n de s�mbolos
NumZeros = 0;             % Cantidad de ceros para sobremuestreo
SwitchChannel = 'ON';     % ON: Enciende el canal de ruido gaussiano,
                          % OFF: Apaga el canal de ruido gaussiano
EsNodB = 2:10;            % SNR en dB [2 3 4 5 6 7 8 9 10]
EsNo = 10.^(EsNodB/10);
SwitchSignal = 'RandomData'; % 'RandomData'
                            % 'ProbeData'
                            % 'RealData'
% K = 5;                    % 11 Coeficientes
K = 6;                      % 13 Coeficientes
%--------------------------------------------------------------------------
% Configuracion de graficos
Graficos = 'OFF';           % Apaga/enciende graficos de datos transmitidos
Constelation = 'ON';
GraficosEqu = 'OFF';          %Canal
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
    % Canal X (613) con distorsion - Adaptivo 
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
        display('Opci�n no valida')
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
    %Constelacion de simbolos recibidos.   
    %ri =  11 ; % Retardo hasta primer simbolo L*T+1 % 1 para adp
    %cSymb = 1000;
    %figure('Name',['Constelacion para ' int2str(EsNodB(k))])
    %stem(RxPAMSignal(ri:T:cSymb),zeros(1,length(ri:T:cSymb)),'b:*')
    
    %----------------------------------------------------------------------
    % Ecualizador del canal X - Adaptivo
     if strcmp(EQU,'ON') 
            delta  = 1/10000;
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
            c_est(1,K+1) = 1; %[0 0 0 0 0 1 0 0 0 0 0]
            c_est_history = zeros(Nent-2*K,2*K+1);
            
            for i=1:1:Nent-2*K 
                y_k = RxPAMSignal(i:i+2*K);  
                z_k(i) = c_est*y_k.';
                e_k = Symbols(i) - z_k(i);      %Error
                c_est_history(i,:) = [c_est];
                c_est = c_est + delta*e_k*y_k;  %Algoritmo LMS
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
                if (mod(k,2)~=0)
                    % Imprimimos respuesta del EQ y respuesta a la salida del EQ
                    figure('Name',['Para ' int2str(EsNodB(k))])
                    subplot(1,2,1) 
                    stem(flip(c_est));        	
                    title('Coeficientes')
                    outEQ = conv(X, [flip(c_est) zeros(1, retardo)]); % Ajustar los ceros segun K
                    subplot(1,2,2) 
                    stem(outEQ(1:1:end));
                    title('Pulso equalizado')
                end
            
                % constelacion de simbolos recibidos.   
%                 ri =  1 ; % Retardo hasta primer simbolo L*T+1 % 1 para adp
%                 cSymb = 200;
%                 figure('Name',['Constelacion para ' int2str(EsNodB(k))])
%                 stem(RxPAMSignal(ri:T:cSymb),zeros(1,length(ri:T:cSymb)),'b:*')                

%             for i=1:2*K+1
%                 subplot(K+1,2,i) 
%                 plot(c_est_history([1:10:100000],i))
%                 coef = i-(K+1);
%                 title(['Coeficiente ' int2str(coef)])
% 
%             end
     end
    
    %----------------------------------------------------------------------
    % Retardo para tener en cuenta los transitorios producidos por las 
    % conv al inicio y al final de la transmisi�n
    RetardoIni = 1; % floor((Sps*Span*2)/2+1) Por 2 de c/filtro /2+1 centrarlo en mx
    ErrorRetardoIni = 0;   % Modificar su valor entre -9 y 9 y justificar degradaci�n respecto a la curva te�rica
    ErrorFrecSimbol = 0;   % Modificar su valor entre 0 y -3 y justificar degradaci�n respecto a la curva te�rica
    Detector = 2*(RxPAMSignal(RetardoIni+ErrorRetardoIni:NumZeros+1+ErrorFrecSimbol:end)>0)-1; % convierte las muestras de los s�mbolos recibidas en valores +-1
    Detector = Detector(1:length(SymbolsAux));
    if(strcmp(Graficos,'ON')&&k==SeleccionGrafico)
        figure
        stem(Detector(1:20));grid;
        title('SENAL DETECTADA');
    end
    NumError(k) = sum(SymbolsAux~=Detector); %Contador de errores
    Rxbits = Detector>0;  %Conversor de s�mbolos a bits
    
    % Constelacion de simbolos recibidos.  
    if (strcmp(Constelation,'ON')&&mod(k,2)~=0)
        ri =  11 ; % Retardo hasta primer simbolo L*T+1 % 1 para adp
        cSymb = 1000;
        figure('Name',['Constelacion para ' int2str(EsNo(k))])
        stem(RxPAMSignal(ri:T:cSymb),zeros(1,length(ri:T:cSymb)),'b:*')
    end
end

%--------------------------------------------------------------------------
%Estimaci�n del tiempo utilizado para la transmisi�n de la informaci�n
M = 2;     %Tipo de modulacion 2PAM, 4PAM, ect
TotalTime = length(Rxbits)/(FreqSymbol*log2(M)*60);
%--------------------------------------------------------------------------
%An�lisis del desempe�o del sistema de comunicaciones mediante curvas de Bit Error Rate (BER):
if strcmp(BERDiagram,'ON')
    simSER = NumError/length(SymbolsAux);  %Estimaci�n del BER obtenido
    figure
    semilogy(EsNodB,0.5*erfc(sqrt(EsNo)),'k'); %Curva teorica
    hold on; semilogy(EsNodB,simSER,'b:o')     %Curva Simulada
    title('Curva de tasa de error de S�mbolos (SER)')
    xlabel('Es/No (dB)')
    ylabel('Probabilidad de error')
    grid on
end