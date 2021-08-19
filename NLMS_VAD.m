clc; clear all;	close all;

%------------SELECIONANDO O ARQUIVO-----------------
% arquivo = 'teste.wav';
 arquivo = 'female_src_1.wav';
% arquivo = 'male_src_1.wav';

[arraySinalRecebido,Fs] = audioread(arquivo);
sizeSinalRecebido=size(arraySinalRecebido);

[arrayFala,Fs] = audioread('male_src_1.wav');

%----------------CRIANDO ECO ARTIFICALMENTE---------
delay = 0.02;
amp = 1;
arrayEco = echo_gen(arraySinalRecebido, Fs, delay, amp);
%---------------------------------------------------

arrayCaptado = arrayEco + [arrayFala;zeros(320,1)];
sizeCaptado = size(arrayCaptado);
 
%--------------------PARAMETROS DO PROGRAMA-------------------

num_ite = 130000; %número de iterações
num_coef = 600;    %número de coeficientes

arrayErro= zeros(1,num_ite); %inicializando array de erro 

for t = 1 : 10          %loop para testar vários parâmetros 
    
    %num_coef = num_coef + 50;
 
	mu = 1; 
	delta = 10^(-4);     %parâmetros do NLMS
	sv = 0.01;             %desvio padrão do ruído aditivo
	MSE = zeros(num_ite,1);      %cria uma matriz de zeros
	arrayFiltrado = zeros(sizeCaptado);

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

	x = arraySinalRecebido;
	z = arrayCaptado;
%	w = rand(M,1);        %cria a matriz coeficientes 
	w = (zeros(1, num_coef + 1))';
    
	%Determinando limiar inicial do VAD.
	Et = z(1:VadLen).^2;
	for(k=1:VadLen)
	    El=El+Et(k);
	end
	%Fim da determinação do limiar inicial do VAD.
    
%-----------------------ITERAÇÕES--------------------------------------------------

	for(i=1:num_ite)        %loop do num de coeficientes até num de iterações 	   
		u=1/norm(x(i+Ns:i+num_coef+Ns))^2/2;
	    %xk = x(i:-1:i-num_coef+1);   %normalizando
	    xk = x(i+Ns:i+num_coef+Ns);
	    y = w'*xk;          %filtro 
	    %e = z(i)- y;          % erro = (sinal filtrado + ruído) - saída do filtro adaptativo 
	    e = z(i) - y;
	    %y=w'*x(i+Ns:i+num_coef+Ns); %Formando a extimativa do ruído (y).
	    %e = z(i+num_coef) - y; %Subtraindo do sinal a extimativa do ruído.

	    %Determinando a energia de uma janela de sinal.
	    Et = z(i:i+num_coef).^2;

        for(j=1:num_coef)
	       E=E+Et(j);
        end
	    %Fim da determinação da energia de uma janela

	    if(i>1)%Quando i é 1, ainda não há o vetor Buffer.
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
%disp('6');
	       	VAD(i)=0;
	     	if(i>1)
	       		if(Var_old~=0)%Condição para não gerar resultado infinito na linha logo abaixo.
	         		Var_relation(i)=var(Buffer)/Var_old;%Calculamos a razão entre a variância do Buffer após adicionarmos a última energia calculada e a variância calculada antes (em ***) e armazenada em Var_old.         
	         	%Dependendo de Var_relation, tomamos uma dentre quatro decisões para atualizar o valor de p, necessário para atualizar o limiar de decisão entre janelas com ou sem voz ativa.
	         	if(Var_relation(i)>=1.25)
	               p=0.25;pes(i)=p;
	           	elseif(Var_relation(i)>=1.10 && Var_relation(i)<=1.25)
	               p=0.20;pes(i)=p;
	           	elseif(Var_relation(i)>=1.00 && Var_relation(i)<=1.10)          
	               p=0.15;pes(i)=p;
	           	elseif(Var_relation(i)<=1.00)
	               p=0.10;pes(i)=p;
	        	end
	       	El=(1-p)*El+p*E;%Atualizando o limiar de decisão entre janelas com ou sem voz ativa.
	    	end
            end        
	    else
	        VAD(i)=1;
	    end
	    Energias(i)=E;
	    E=0;%Voltamos E para 0 para poder reiniciar os cálculos de energia.
	    Limiares(i)=El;    		

	MSE(i) = MSE(i)+e^2;
		    
	arrayErro(i) = e;

	end
	arrayErro = arrayErro'; %transpondo array erro 

% ax1 = nexttile;
% plot(ax1,arrayErro(:,1))
% title(ax1,'Filtrado - Erro')

%figure();
%plot(w);

% ax2 = nexttile;
% plot(ax2,10*log10(MSE),'r')
% title(ax2,'r')



end

