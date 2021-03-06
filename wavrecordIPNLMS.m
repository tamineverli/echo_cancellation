clc;%Limpa a tela.
clear all;%Apaga todas as variáveis existentes no workspace.
close all;%Fecha todas as janelas abertas pelo matlab.

arquivo = 'male_src_1.wav';
[arrayEntrada,Fs] = audioread(arquivo);
sizeArrayEntrada=size(arrayEntrada);
arrayEntrada = pinknoise(sizeArrayEntrada);
%arrayEntrada = randn(sizeArrayEntrada);


%RIR do código do Matheus a partir daqui 
%ajuste de parâmetros da RIR 
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

%--------------------PARAMETROS DO PROGRAMA-------------------

N = 160000; %número de iterações
M = 300;   %número de coeficientes
p=0;

    Nd = 30; % Tamanho do canal a ser identificado 
    %N=length(arrayCaptado)-100; %número de iterações

    beta   = 0.20;      % beta (step size paramenter)
    ro     = 0.01;      % ro
    delta  = 0.01;      % delta
    deltap = 0.01;      % deltap
    alpha  = -0.5;

    L=1;    %número de realizações
    mu=1; 
    delta=10^(-4);  %parâmetros do NLMS
    sv=0.0001;  %desvio padrão do ruído aditivo
    MSE=zeros(N,1);      %cria uma matriz de zeros
    MSD = zeros(N,1);
    arrayFiltrado=zeros(sizeCaptado);

    inactiv_it = 600; % Iterations without algorithm activated

    VadLen=160;%Tamanho da janelas de VAD para calculo de limiar inicial.
    E=0;  %Energia da janela.
    El=0; %Energia limiar do VAD.
    LMult=1.5;%Fator de multiplicação do limiar. Para baixo SNR, deixar próximo de 1.0
    Et=zeros(1,VadLen); %Variável de cálculo temporário.
    Buffer_pos=1;%Inicia o índice que armazenará as energias no vetor Buffer.
    BuffLen=40;%Tamanho do Buffer de VAD.
    Var_relation=0;%Inicia a variável `relação de variâncias`, que será usada pelo VAD mais tarde.
    Var_old=0;%Inicia a variável `variância old`, que será usada pelo VAD mais tarde.
    p=0.2;%Inicia a variável p, que será utilizada para atualizar o limiar de detecção do VAD.

%--------------------REALIZAÇÕES--------------------

    for k=1:L            %loop no número de realizações 
       x = arrayEntrada; 
       d = arrayCaptado;
       wg = rand(M,1); 
            
%-----------------------ITERAÇÕES--------------------------------------------------
        
            for i = inactiv_it : N
                % CALCULO DO ERRO
                xi = x(i:-1:i-M+1); %% TRANSPOR SE PRECISAR 
                yg = wg' * xi;    
                eg = d(i) - yg;
          
%----------------------ADAPTAÇÃO---------------------------------------------------

                 if i > inactiv_it
                    Fw = Fseg( abs(wg) ); % mi-law
                    Fw = ( abs(wg) )'; % sem mi-law
                    sumh=sum(Fw); 
                    gama_min = ro * deltap;
                    gden = 0;   % denominador de g
                   

                    for n = 1 : M
                        gama( n ) = (1-alpha)*sumh/M+(1+alpha)*Fw(n);
                        % montando parte do denominador de g
                        gden = gden + gama( n );
                        gama(n)=max([gama(n) gama_min]);
                       
                
                    end

                    gden=max([gden delta]);
                    g = gama ./ gden;
                    G = diag( g );
                    auxg = x(i:-1:i-M+1); %% TRANSPOR SE PRECISAR 
                    
                    numerador = (beta * G * auxg * eg);
                    denominador = (auxg' * G * auxg + delta );
                    wg = wg +  numerador / denominador ;  % new coefficients of the estimated filter             
            end %if i > inactiv_it
         arrayFiltrado(i)=eg; 
         MSE(i) = MSE(i)+eg^2;
         MSD(i) = sum((wg' - h_air(1,1:M)).^2); 
         end
   
    end

figure
plot(10*log10(MSE))
title('MSE - White Noise')
figure
plot(10*log10(MSD))
title('MSD - White Noise')
