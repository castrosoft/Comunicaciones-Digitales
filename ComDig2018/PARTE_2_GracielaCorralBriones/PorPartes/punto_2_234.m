% Configurar el simulador para estudiar el efecto de la interferencia 
% intersímbolo. Para ello reemplazar el los módulos conformador, canal con 
% ruido gaussiano y filtro apareado por el canal del ejemplo 6.11 del Libro
% de Proakis. Nota: ajustar los tasa de sobremuestreo para que sea 
% compatible con la respuesta del canal
% Analizar el efecto de la interferencia intersímbolo (anular el ruido) en 
% la constelación de símbolos recibidos. 
clear all;close all;clc;
%Parámetros
DataLength = 1000;        % Cantidad de símbolos a transmitir
FreqSymbol = 1.2e3;         % Frecuencia de transmisión de símbolos
NumZeros = 1;               % Cantidad de ceros para sobremuestreo
SwitchChannel = 'ON';      % ON: Enciende el canal de ruido gaussiano,
                            % OFF: Apaga el canal de ruido gaussiano
EsNodB = 2:9;               % SNR en dB
EsNo = 10.^(EsNodB/10);
SwitchSignal = 'RandomData';% 'RandomData'
                            % 'ProbeData'
                            % 'RealData'
K = 6;                      % 11 Coeficientes
% K = 6;                      % 13 Coeficientes
%--------------------------------------------------------------------------
% Configuracion de graficos
Graficos = 'OFF';           % Apaga/enciende graficos de datos transmitidos
Constelation = 'ON';
GraficosEqu = 'ON';
BERDiagram = 'ON';
SeleccionGrafico = 8;       % Determina grafico segun k
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
    % Canal H con distorcion 

    %   Ajustes del canal
    T  = NumZeros + 1;
    Ts = T/2; 
    L = 5;
    t = -L*T:Ts:L*T;
    H = 1./(1+((T/2)*t).^2);

    if(strcmp(GraficosEqu,'ON')&&k==SeleccionGrafico)
        figure
        stem(H);
        title('RESPUESTA AL IMPULSO DEL CANAL H');
    end
    
    % Se pasa simbolos por filtro transmisor
    TxPAMSignal = conv(Symbols,H);
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
    % Ecualizador del canal H con ZF
    
    % Matriz X
    MX = zeros((2*K+1),(2*K+1));
    for m = -K:1:K
        for n = -K:1:K
            MX(m+K+1,n+K+1) = (1/(1+(2*(m*T-n*Ts)/T)^2));
        end
    end
    qSol  = zeros(2*K+1,1);
    qSol(K+1,1) = 1;
    c_opt = inv(MX)*qSol; %c optimos del ZF
    display(c_opt);
    % Imprimimos respuesta del EQU y respuesta a la salida del EQU
    if(strcmp(GraficosEqu,'ON')&&k==SeleccionGrafico)
        figure
        stem(c_opt);
        title('RESPUESTA DEL EQU ZF(Coeficientes)');
        figure
        outEQ = conv(c_opt,H);
        % Downsampling
        if K==5
            stem(outEQ(2:2:end)); %Para K=5
        else
            stem(outEQ(1:2:end)); %Para K=6
        end
        title('RESPUESTA A LA SALIDA DEL EQU');
    end
    outFilter = conv(c_opt,RxPAMSignal);
    RxPAMSignal = outFilter(K+1:end);
    
    %----------------------------------------------------------------------
    % Retardo para tener en cuenta los transitorios producidos por las 
    % conv al inicio y al final de la transmisión
    RetardoIni = 11; % floor((Sps*Span*2)/2+1) Por 2 de c/filtro /2+1 centrarlo en mx
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
    if (strcmp(Constelation,'ON')&&k==SeleccionGrafico)
        ri =  11 ; % Retardo hasta primer simbolo L*T+1 % 1 para adp
        cSymb = 1000;
        figure('Name',['Constelacion para ' int2str(EsNo(k))])
        stem(RxPAMSignal(ri:T:cSymb),zeros(1,length(ri:T:cSymb)),'b:*')
    end
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