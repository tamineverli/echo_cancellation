clear all

arquivo = 'female_src_1.wav';
arquivoCaptado = 'femaleCaptado.wav';

% arquivo = 'male_src_1.wav';
% arquivoCaptado = 'maleCaptado.wav';

[arrayEntrada,Fs] = audioread(arquivo);
sizeArrayEntrada=size(arrayEntrada);

%[arrayCaptado,Fs] = audioread(arquivoCaptado);
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
[h_air,air_info] = LoadAIR.loadAIR(airpar,'AIR_LIB\');


load rir.mat
 size(h_air)
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
        e = d(k)- y;      % erro = (sinal filtrado + ruído) - saída do filtro adaptativo 
        w = w + mu*e*xk /((xk'*xk)+delta); %atualizacao dos coeficientes 
        matrizCoeficientes(:,k) = w;
        MSE(k) = MSE(k)+e^2;
        arrayY (k) = y;
        arrayErro = e;
        arrayFiltrado(k) = e;
        
        
    end
    
    
norm_h_air = norm(h_air).^2;

for i = 1:N
    norm_wh = norm(matrizCoeficientes(:,i) - h_air).^2;
    MSD (i) = (norm_wh)/norm_h_air; 
end


 MSE = (MSE - min(MSE))/(max(MSE) - min(MSE)); 

%  plot(10*log10(MSD),'b')

%title(ax5,'r')
%figure
plot(10*log10(MSE),'r')
%title(ax5,'r')
