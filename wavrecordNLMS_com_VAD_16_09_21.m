clear all
close all

arquivo = 'male_src_1.wav';
[arrayEntrada,Fs] = audioread(arquivo);
sizeArrayEntrada=size(arrayEntrada);
%arrayEntrada = pinknoise(sizeArrayEntrada);
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
[h_air,air_info] = LoadAIR.loadAIR(airpar,'AIR_LIB\');

%load rir.mat
%size(h_air)
%h_air = h_air/norm(h_air);
%figure
%plot(h_air)
%RIR=randn(1000,1);
arrayCaptado=conv(arrayEntrada,h_air);
sizeCaptado=size(arrayCaptado);


N = 160000; %Numero de Iteracoes 
M = 4000;   %Numero de coeficientes

mu=1; 
delta = 10^(-4);     %parÃ¢metros do NLMS
sv = 0.01;             %desvio padrÃ£o do ruÃ­do aditivo
MSE = zeros(N,1);      %cria uma matriz de zeros
MSD = zeros(N,1);

%----------------Parametros VAD------------------------

VadLen = 160;		%Tamanho da janelas de VAD para calculo de limiar inicial.
E = 0;  			%Energia da janela.
El = 0;				%Energia limiar do VAD.
%LMult = 1.5        %Fator de multiplicação do limiar. Para baixo SNR, deixar próximo de 1.0
LMult = 20;

Et = zeros(1,VadLen); %Variável de cálculo temporário.
Buffer_pos = 1;		%Inicia o índice que armazenará as energias no vetor Buffer.
BuffLen = 40;		%Tamanho do Buffer de VAD.
Var_relation = 0; 	%Inicia a variável `relação de variâncias`, que será usada pelo VAD mais tarde.
Var_old = 0;		%Inicia a variável `variância old`, que será usada pelo VAD mais tarde.
p = 0.2;			%Inicia a variável p, que será utilizada para atualizar o limiar de detecção do VAD.
Buffer = zeros(1,BuffLen);

Ns=142;%Ajuste da correlação dos sinais
  
%-------------------------------------------------------------
arrayFiltrado = zeros(sizeCaptado);

x = arrayEntrada;
d = arrayCaptado;
w = zeros(M,1);        %cria a matriz coeficientes 

%-----------Determinando limiar inicial do VAD----------------
	
Et = d(1:VadLen).^2;
	for(k=1:VadLen)
	    El=El+Et(k);
	end
%-------------------------------------------------------------
  

    for k = M:N         
       
        xk=x(k:-1:k-M+1);   %normalizando 
        y = w'*xk;
        e = d(k)- y;      % erro = (sinal filtrado + ruido) - saida do filtro adaptativo 
     %  w = w + mu*e*xk /((xk'*xk)+delta); %atualizacao dos coeficientes 
        
        %Determinando a energia de uma janela de sinal.
	      Et = d(k:-1:k-M+1).^2; %Preciso confirmar com Mariane 

        for(j=1:M)
	       E=E+Et(j);
        end
	    %Fim da determinação da energia de uma janela

	    if(k>M)%Quando i é 1, ainda não há o vetor Buffer.
	    	Var_old=var(Buffer);%(***)Determina a variância do Buffer antes de inserir a última energia calculada.
	    end

	    Buffer(Buffer_pos)=E;%Armazena a energia calculada como último elemento do Buffer.
	    Buffer_pos=Buffer_pos+1;%Atualiza o índice do Buffer para a próxima iteração.
	    
	    if(Buffer_pos==BuffLen)%Forçamos o Buffer a ter BuffLen elementos.
	        Buffer_pos=1;
	    end


	    if(E<LMult*El) %Verifica, pela potência do sinal, se há voz. Havendo, não atualizamos os coeficientes.
	       	%w=w+u*err(k)*z(i+Ns:i+N+Ns);%Atualiza os coeficientes do filtro.
	       	w = w + mu*e*xk /((xk'*xk)+delta); %atualização dos coeficientes 

	       	VAD(k)=0;
	     	if(k>M)
	       		if(Var_old~=0)%Condição para não gerar resultado infinito na linha logo abaixo.
	         		Var_relation(k)=var(Buffer)/Var_old;%Calculamos a razão entre a variância do Buffer após adicionarmos a última energia calculada e a variância calculada antes (em ***) e armazenada em Var_old.         
	         	%Dependendo de Var_relation, tomamos uma dentre quatro decisões para atualizar o valor de p, necessário para atualizar o limiar de decisão entre janelas com ou sem voz ativa.
	         	if(Var_relation(k)>=1.25)
	               p=0.25;pes(k)=p;
	           	elseif(Var_relation(k)>=1.10 && Var_relation(k)<=1.25)
	               p=0.20;pes(k)=p;
	           	elseif(Var_relation(k)>=1.00 && Var_relation(k)<=1.10)          
	               p=0.15;pes(k)=p;
	           	elseif(Var_relation(k)<=1.00)
	               p=0.10;pes(k)=p;
	        	end
	       	El=(1-p)*El+p*E;%Atualizando o limiar de decisão entre janelas com ou sem voz ativa.
	    	end
        
     end
     
        
          
	    else
	        VAD(k)=1;
	    end
      Energias(k)=E;
	    E=0;%Voltamos E para 0 para poder reiniciar os cálculos de energia.
	    Limiares(k)=El;
      
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