clear all
close all

arquivo = 'male_src_1.wav';
[arrayEntrada,Fs] = audioread(arquivo);
sizeArrayEntrada=size(arrayEntrada);
%arrayEntrada = pinknoise(sizeArrayEntrada);
%arrayEntrada = randn(sizeArrayEntrada);


%RIR do c�digo do Matheus a partir daqui 
%ajuste de par�metros da RIR 
airpar.fs = 48e3;  
% ## airpar.rir_type = 1;
% ## airpar.room = 4;
% ## airpar.channel = 1;
% ## airpar.head = 1;
% ## airpar.rir_no = 4;
airpar.rir_type = 1;
airpar.room = 2;
airpar.channel = 1;
airpar.head = 0;
airpar.rir_no = 1;
%[h_air,air_info] = LoadAIR.loadAIR(airpar,'AIR_LIB\');

load rir.mat
size(h_air)
h_air = h_air/norm(h_air);
%figure
%plot(h_air)
%RIR=randn(1000,1);
arrayCaptado=conv(arrayEntrada,h_air);
sizeCaptado=size(arrayCaptado);

arrayFiltrado1 = zeros(1,200000);
arrayFiltrado2 = zeros(1,200000);

N = 160000; %número de iterações
M = 4000;   %número de coeficientes

mu=1; 
delta = 10^(-4);     %parâmetros do NLMS
sv = 0.01;             %desvio padrão do ruído aditivo
MSE = zeros(N,1);      %cria uma matriz de zeros
MSD = zeros(N,1);
MSE_Medio = zeros(N,1);
arrayFiltrado = zeros(sizeCaptado);

x = arrayEntrada;
d = arrayCaptado;
w = zeros(M,1);        %cria a matriz coeficientes 
matrizCoeficientes = zeros(M,N);

    for k = M:N        %loop do num de coeficientes até num de iterações 
       
        xk=x(k:-1:k-M+1);   %normalizando 
        y = w'*xk;
        e = d(k)- y;      % erro = (sinal filtrado + ruido) - saida do filtro adaptativo 
        w = w + mu*e*xk /((xk'*xk)+delta); %atualizacao dos coeficientes 
        MSE(k) = MSE(k)+e^2;
        MSD(k) = sum((w' - h_air(1,1:M)).^2); 
        arrayFiltrado(k) = e;
    end    

%MSE = (MSE - min(MSE))/(max(MSE) - min(MSE)); 

%chplot(10*log10(MSD),'b')

figure
plot(10*log10(MSE))
title('MSE - Male')
figure
plot(10*log10(MSD))
title('MSD - Male')
