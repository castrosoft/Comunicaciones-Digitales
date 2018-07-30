%=======================================================================
% Comunicaciones Digitales
% Prof.: Dr. Mario Hueda
% Practico Lab. 2
%=======================================================================
clear
close all
%	2.1 Lectura de las muestras de ruido.
N=50000;
n=1/sqrt(2)*(randn(1,N)+j*randn(1,N));
%	P2.2	Generacion de simbolos binarios.
btxd=randint(1,N);
%btxd=.0*ones(1,N);
x=(btxd-.5)/.5;
hist(x,100);
title('Histograma de los simbolos transmitidos');
%pause;
%	P2.3	SNR de interes = 6.8 dB (por ejemplo, para BER=1e-3).
SNR=6.8;
%SNR=4.3;
SNRv=10^(SNR/10);
No=1/SNRv;
n=sqrt(No)*n;
%	P2.4	Muestra recibida.
r=real(x+n);
%	P2.5	Histograma de las muestras recibidas.
hist(r,100)
title('Histograma de las muestras recibidas');
disp('Pulse una tecla para continuar');
%pause;
%	P2.6 Detector y calculo de errores
i=find(r>0);
brxd=zeros(1,N);
brxd(i)=1;
fprintf('\nError de simulacion = %f',mean(xor(btxd,brxd)));
fprintf('\nError de teoria = %f\n',0.5*erfc(sqrt(SNRv)));

%Nota: Q(x)=0.5*erfc(x/sqrt(2));