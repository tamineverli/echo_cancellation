clear all

arquivo = 'female_src_1.wav';
arquivoCaptado = 'femaleCaptado.wav';

% arquivo = 'male_src_1.wav';
% arquivoCaptado = 'maleCaptado.wav';

[arrayEntrada,Fs] = audioread(arquivo);
sizeArrayEntrada=size(arrayEntrada);

[arrayCaptado,Fs] = audioread(arquivoCaptado);
sizeCaptado=size(arrayCaptado);


arrayFiltrado1 = zeros(1,200000);
arrayFiltrado2 = zeros(1,200000);

N = 160000; %número de iterações
M = 1000;   %número de coeficientes


mu=1; 
delta = 10^(-4);     %parâmetros do NLMS
sv=0.01;             %desvio padrão do ruído aditivo
MSE=zeros(N,1);      %cria uma matriz de zeros
arrayFiltrado=zeros(sizeCaptado);

x = arrayEntrada;
d = arrayCaptado;
w = zeros(M,1);        %cria a matriz coeficientes 

for k = M:N        %loop do num de coeficientes até num de iterações 
   
    xk=x(k:-1:k-M+1);   %normalizando 
    y = w'*xk;
    e = d(k)- y;      % erro = (sinal filtrado + ruído) - saída do filtro adaptativo 
    w = w + mu*e*xk /((xk'*xk)+delta); %atualização dos coeficientes 
    MSE(k) = MSE(k)+e^2;
    
    
    arratY (k) = y;
    arrayErro = e;
    arrayFiltrado(k) = e;
    
    
end


arrayFiltrado = arrayFiltrado';

ax5 = nexttile;
plot(ax5,10*log10(MSE),'r')
title(ax5,'r')
