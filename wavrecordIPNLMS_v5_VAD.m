clc;%Limpa a tela.
clear all;%Apaga todas as variáveis existentes no workspace.
close all;%Fecha todas as janelas abertas pelo matlab.

% arquivo = 'teste.wav';
% arquivoCaptado = 'testeCaptado.wav';
% 
arquivo = 'female_src_1.wav';
arquivoCaptado = 'femaleCaptado.wav';
% 
% arquivo = 'male_src_1.wav';
% arquivoCaptado = 'maleCaptado.wav';

[arrayMusica,Fs] = audioread(arquivo);
sizeArrayMusica=size(arrayMusica);

[arrayCaptado,Fs] = audioread(arquivoCaptado);
sizeCaptado=size(arrayCaptado);


%--------------------PARAMETROS DO PROGRAMA-------------------

N = 160000; %número de iterações
M = 300;   %número de coeficientes
p=0;

for p = 1 : 3
        
    %N = N + 100000;
    M = M + 30;
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
    %h=rand(M,1);         %cria uma matriz rand de coeficientes
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
       x = arrayMusica; 
       d = arrayCaptado;
       wg = rand(M,1); 
            
%-----------------------ITERAÇÕES--------------------------------------------------
        
            for i = inactiv_it : N
                % CALCULO DO ERRO
                xi = x( i:-1:i-M+1); %% TRANSPOR SE PRECISAR 
                yg( i ) = wg' * xi;    
                eg(i)  = d(i)   - yg(i);
          
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
                    
                    numerador = (beta * G * auxg * eg(i));
                    denominador = (auxg' * G * auxg + delta );
                    wg = wg +  numerador / denominador ;  % new coefficients of the estimated filter             
            end %if i > inactiv_it
         end
         arrayFiltrado=eg; 
         MSE=MSE+eg'.^2;
    end

    ax1 = nexttile;
    plot(ax1,arrayCaptado(:,1))
    title(ax1,'Captado')

    ax2 = nexttile;
    plot(ax2, arrayMusica(:,1))
    title(ax2,'Musica')

    ax3 = nexttile;
    plot(ax3,arrayFiltrado(1,:))
    title(ax3,'Filtrado')

    ax4 = nexttile;
    plot(ax4,10*log10(MSE/L),'r')
    title(ax4,'r')
    %disp('FIM')
end
