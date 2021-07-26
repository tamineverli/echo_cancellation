clear all

% arquivo = 'teste.wav';
% arquivoCaptado = 'testeCaptado.wav';

arquivo = 'female_src_1.wav';
arquivoCaptado = 'femaleCaptado.wav';

% arquivo = 'male_src_1.wav';
% arquivoCaptado = 'maleCaptado.wav';

[arrayMusica,Fs] = audioread(arquivo);
sizeArrayMusica=size(arrayMusica);

[arrayCaptado,Fs] = audioread(arquivoCaptado);
sizeCaptado=size(arrayCaptado);

% clear all
% [arrayMusica,Fs] = audioread('male_src_1.wav');
% sound(arrayMusica,Fs);
% sizeArrayMusica=size(arrayMusica);
% objCaptado= audiorecorder(Fs,8,1);
% disp('Start Recording.')
% recordblocking(objCaptado, 10);
% disp('End of Recording.');
% arrayCaptado = getaudiodata(objCaptado);
% sound(arrayCaptado,Fs);
% sizeCaptado=size(arrayCaptado);
% filename = 'maleCaptado.wav';
% audiowrite(filename,arrayCaptado,Fs);

arrayFiltrado1 = zeros(1,200000);
arrayFiltrado2 = zeros(1,200000);

N = 160000; %número de iterações
M = 300;   %número de coeficientes
p=0;

for p = 1 : 10
    
    %N = N + 100000;
    M = M + 30;
 
mu=1; 
delta = 10^(-4);     %parâmetros do NLMS
sv=0.01;             %desvio padrão do ruído aditivo
MSE=zeros(N,1);      %cria uma matriz de zeros
arrayFiltrado=zeros(sizeCaptado);

x = arrayMusica;
d = arrayCaptado;
w=zeros(M,1);        %cria a matriz coeficientes 

%ir aumentando M, começar por uns 300 por conta da reverberação 
%ir aumentando o número de interações, na casa dos 10k

for k=M:N        %loop do num de coeficientes até num de iterações 
   
    xk=x(k:-1:k-M+1);   %normalizando 
    y = w'*xk;
    e=d(k)- y;      % erro = (sinal filtrado + ruído) - saída do filtro adaptativo 
    w= w + mu*e*xk /((xk'*xk)+delta); %atualização dos coeficientes 
    MSE(k)=MSE(k)+e^2;
    
    
    arratY (k) = y;
    arrayErro = e;arrayFiltrado1(k) = e;
    arrayFiltrado2(k) = d(k);
    
end


arrayFiltrado1 = arrayFiltrado1';
arrayFiltrado2 = arrayFiltrado2';
 
% ax1 = nexttile;
% plot(ax1,arrayCaptado(:,1))
% title(ax1,'Captado')
% 
% ax2 = nexttile;
% plot(ax2, arrayMusica(:,1))
% title(ax2,'Musica')

ax3 = nexttile;
plot(ax3,arrayFiltrado1(:,1))
title(ax3,'Filtrado - Erro')

% ax4 = nexttile;
% plot(ax4,arrayFiltrado2(:,1))
% title(ax4,'Filtrado - "D"')

% ax5 = nexttile;
% plot(ax5,10*log10(MSE),'r')
% title(ax5,'r')

%sound(arrayFiltrado1, Fs);
end