clear all

% arquivo = 'teste.wav';
% arquivoCaptado = 'testeCaptado.wav';
% % 
arquivo = 'female_src_1.wav';
arquivoCaptado = 'femaleCaptado.wav';
% % 
% % arquivo = 'male_src_1.wav';
% % arquivoCaptado = 'maleCaptado.wav';

[arrayMusica,Fs] = audioread(arquivo);
sizeArrayMusica=size(arrayMusica);

[arrayCaptado,Fs] = audioread(arquivoCaptado);
sizeCaptado=size(arrayCaptado);


%--------------------------------------------
% PARAMETROS DO PROGRAMA
%-------------------------------------------

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
%h=rand(M,1);         %cria uma matriz rand de coeficientes
arrayFiltrado=zeros(sizeCaptado);

% Iterations without algorithm activated
inactiv_it = 600;

%--------------------------------------------
%REALIZAÇÕES
%disp('REALIZAÇÕES')
%-------------------------------------------

for k=1:L            %loop no número de realizações 
   x = arrayMusica; 
   d = arrayCaptado;
   wg = rand(M,1); 
        
%-------------------------------------------------------------------------
% ITERAÇÕES
%disp('ITERAÇÕES')
%-------------------------------------------------------------------------

        for i = inactiv_it : N
            % CALCULO DO ERRO
            xi = x( i:-1:i-M+1); %% TRANSPOR SE PRECISAR 
            yg( i ) = wg' * xi;    
            eg(i)  = d(i)   - yg(i);
      
%-------------------------------------------------------------------------
% ADAPTAÇÃO
%disp('ADAPTAÇÃO')
%-------------------------------------------------------------------------        
             if i > inactiv_it
                Fw = Fseg( abs(wg) ); % mi-law
                Fw = ( abs(wg) )'; % sem mi-law
                sumh=sum(Fw); 
                gama_min = ro * deltap;
                gden = 0;   % denominador de g
               
%disp('ADAPTAÇÃO1')
                for n = 1 : M
                    gama( n ) = (1-alpha)*sumh/M+(1+alpha)*Fw(n);
                    % montando parte do denominador de g
                    gden = gden + gama( n );
                    gama(n)=max([gama(n) gama_min]);
                   
%disp('ADAPTAÇÃO2')                
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

ax4 = nexttile;
plot(ax4,10*log10(MSE/L),'r')
title(ax4,'r')
%disp('FIM')
