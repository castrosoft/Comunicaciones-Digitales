% Configurar el simulador para estudiar el efecto de la interferencia 
% intersímbolo. Para ello reemplazar el los módulos conformador, canal con 
% ruido gaussiano y filtro apareado por el canal del ejemplo 6.11 del Libro
% de Proakis. Nota: ajustar los tasa de sobremuestreo para que sea 
% compatible con la respuesta del canal
% Analizar el efecto de la interferencia intersímbolo (anular el ruido) en 
% la constelación de símbolos recibidos. 
clear all;close all;clc;
%Parámetros
DataLength = 500;        % Cantidad de símbolos a transmitir
FreqSymbol = 1.2e3;         % Frecuencia de transmisión de símbolos
NumZeros = 1;               % Cantidad de ceros para sobremuestreo
SwitchChannel = 'OFF';      % ON: Enciende el canal de ruido gaussiano,
                            % OFF: Apaga el canal de ruido gaussiano
EsNodB = 2:9;               % SNR en dB
EsNo = 10.^(EsNodB/10);
SwitchSignal = 'RandomData';% 'RandomData'
                            % 'ProbeData'
                            % 'RealData'
%--------------------------------------------------------------------------
% Configuracion de graficos
Graficos = 'ON';           % Apaga/enciende graficos de datos transmitidos
Constelation = 'ON';
SeleccionGrafico = 1;       % Determina grafico segun k
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

    if(strcmp(Graficos,'ON')&&k==SeleccionGrafico)
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