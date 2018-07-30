%1) Configurar el simulador para cumplir criterio de Nyquist.
%1.1 Graficar diagrama ojo a la salida del filtro apareado para distintos roll offs
clear all;close all;clc;
%Parámetros
DataLength = 1000000;       % Cantidad de símbolos a transmitir
FreqSymbol = 1.2e3;         % Frecuencia de transmisión de símbolos
NumZeros = 9;              % Cantidad de ceros para sobremuestreo
SwitchChannel = 'ON';       % ON: Enciende el canal de ruido gaussiano,
                            % OFF: Apaga el canal de ruido gaussiano
EsNodB = 2:10;               % SNR en dB
EsNo = 10.^(EsNodB/10);
SwitchSignal = 'RandomData';% 'RandomData'
                            % 'ProbeData'
                            % 'RealData'
RxFilter = 'Matched';       %'Matched'
                            %'NotMatched'
%--------------------------------------------------------------------------
% Diseño del filtro conformador de la señal
RollOff = 0.5;
Span = 10;  % Esto es el span del raised cos en duracion de symbolo
Sps = 10;   % Muestras por span, o sea, cada 10 (Sps) debe ir un simnbolo

SRRCFilter = rcosdesign(RollOff,Span,Sps,'sqrt'); 
RCFilter = rcosdesign(RollOff,Span,Sps,'normal'); 
stem(SRRCFilter);
%--------------------------------------------------------------------------
% Configuracion de graficos
Graficos = 'OFF';           % Apaga/enciende graficos de datos transmitidos
EyeDiagram = 'OFF';
BERDiagram = 'ON';
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
    
    % Se pasa simbolos por filtro transmisor
    TxPAMSignal = conv(Symbols,SRRCFilter);
    %----------------------------------------------------------------------
    %Canal Gausiano:
    if strcmp(SwitchChannel,'ON')
        Noise = randn(1,length(TxPAMSignal))./sqrt(2*EsNo(k));
    elseif strcmp(SwitchChannel,'OFF')
        Noise = zeros(1,length(TxPAMSignal));
    else
        display('Opción no valida')
    end
    NoisySignal = TxPAMSignal+Noise;
    %----------------------------------------------------------------------
    %Receptor Digital:
    if strcmp(RxFilter,'Matched')
        RxPAMSignal = conv(NoisySignal,SRRCFilter);
    elseif strcmp(RxFilter,'NotMatched')
        RxPAMSignal = conv(NoisySignal,RCFilter);
    end
    if(strcmp(Graficos,'ON')&&k==SeleccionGrafico)
        figure
        stem(RxPAMSignal(1:200));grid;
        title('SENAL RECIBIDA DESPUES DEL FILTRO, CON SOBREMUESTREO');
    end

    % Retardo para tener en cuenta los transitorios producidos por las 
    % conv al inicio y al final de la transmisión
    RetardoIni = 101; % floor((Sps*Span*2)/2+1) Por 2 de c/filtro /2+1 centrarlo en mx
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
end

%--------------------------------------------------------------------------
% Dibuja diagrama de ojo
if strcmp(EyeDiagram,'ON')
    eyediagram(RxPAMSignal(RetardoIni:2000+RetardoIni),NumZeros+1);grid on;
    title(['Diagrama de Ojo para RollOff: ' num2str(RollOff)]);
    pause(2)
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