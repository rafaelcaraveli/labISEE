%========        PRATICA 7 ISEE - 2021     ========
%========  Atividade 1 - Fluxo de Potência ========
%========            Abril/2021            ========
%========  Amanda; Bruno; Thiago ; Rafael  ========

% Descrição: Sistesma com 6 barras sendo:
% Barra 1: Referência
% Barra 2: PV (geração)
% Barra 3: PV (geração)
% Barra 4: PQ (carga)
% Barra 5: PQ (carga)
% Barra 6: PQ (carga)

clear;
clc;
close all;

%% Valores Base
Sb=100e6;

%%  ------- Informações Conhecidas  -------
V=[1.02 1.04 1.01 0 0 0];  %Connhece V de 1,2,3
theta=[0 0 0 0 0 0];       %Só é conhecido o angulo de 1

% Para potência ativa foi utilizado os valores referentes a letra a
Pk=[0 1.5 1.5 -1 -1 -1];  %Conhece P de 2 a 6


% Potência reativa
fp=[0.9806 0.8944 0.9950];
t_aux=tan(acos(fp)).*Pk(4:6);
Qk=[0 0 0 t_aux(1) t_aux(2) t_aux(3)]; %Conhece 4 a 6


% Montando Ybus 
z=0.04+0.08j; %[pu]
%z=0.00+0.08j; %Sem perdas
y=1/z;
ysh=0.02*j; %Susceptância shunt da LT




Ybus=[2*y+2*ysh/2, -y, 0, 0, -y, 0;
      -y,  2*y+2*ysh/2,0, -y, 0, 0;
      0, 0, 2*y+2*ysh/2, 0, -y, -y;
      0, -y, 0, 3*y+3*ysh/2, -y, -y;
      -y, 0, -y, -y, 3*y+3*ysh/2, 0;
      0, 0, -y, -y, 0, 2*y+2*ysh/2];
G=real(Ybus);
B=imag(Ybus);


%% ------- Estimativa de valores iniciais -------
%Para as tensões não conhecidas (4 a 6) estipulou-se o valor inicial de 1pu
%Para os ang não conhecidos (2 a 6) utilizou os valores do fluxo CC (a)
V=[1.02 1.04 1.01 1 1 1];
theta=[0 0.03733 0.0133 -0.04533 -0.03733 -0.05600]; 

%% ------- Cálculos para o subsistema 1  -------
% Nessa etapa é cálculdo Vk thetak para as barras PQ e thetak para as barras PV
% Se calcula deltaP para as barras PV e PQ (2 a 6) e deltaQ para as barras PQ (4 a 6)
% Os valores já são inseridos diretamente em g(x)=[deltaP;deltaQ] 

dP2=Pk(2)-V(2)*sum(V.*(G(2,:).*cos(theta(2)-theta) + B(2,:).*sin(theta(2)-theta)));
dP3=Pk(3)-V(3)*sum(V.*(G(3,:).*cos(theta(3)-theta) + B(3,:).*sin(theta(3)-theta)));
dP4=Pk(4)-V(4)*sum(V.*(G(4,:).*cos(theta(4)-theta) + B(4,:).*sin(theta(4)-theta)));
dP5=Pk(5)-V(5)*sum(V.*(G(5,:).*cos(theta(5)-theta) + B(5,:).*sin(theta(5)-theta)));
dP6=Pk(6)-V(6)*sum(V.*(G(6,:).*cos(theta(6)-theta) + B(6,:).*sin(theta(6)-theta)));
dQ4=Qk(4)-V(4)*sum(V.*(G(4,:).*sin(theta(4)-theta) - B(4,:).*cos(theta(4)-theta)));
dQ5=Qk(5)-V(5)*sum(V.*(G(5,:).*sin(theta(5)-theta) - B(5,:).*cos(theta(5)-theta)));
dQ6=Qk(6)-V(6)*sum(V.*(G(6,:).*sin(theta(6)-theta) - B(6,:).*cos(theta(6)-theta)));

gx=[dP2; dP3; dP4; dP5; dP6; dQ4; dQ5; dQ6];
%------------

% Para obter o resultado é utilizado o método de Newton-Raphson, o objetivo
% é que a diferença entre a potência que já conhecemos e o somatorio das
% potências saindo da barra seja 0. Logo, sai do loop quando o erro do 
% maior termo de gx (em módulo) for menor ou igual que o erro estipulado.
% Se realizar 20 interações e não obter o valor considera-se que o sistema
% não convergiu.

%Erro Estipulado
erro=0.003e-5;
%erro=0.003;
iteracao=20; 

for it=0:iteracao
    itp(it+1)=it+1       %variavel para plotar
    errop(it+1)=max(abs(gx)); %variavel para plotar
   % fprintf('%d - %f\n',it,max(abs(gx)));
   
    if max(abs(gx)) > erro
    % Jc = [H N
    %       M L];
    % Para H -> temos P em relação a theta => Muito sensível
    % Para N -> temos P em relação a tensão => Pouco sensível
    % Para M -> temos Q em relação a theta => Pouco sensível
    % Para L -> temos Q em relação a tensão => Muito sensível
    for k = 2:6                                         % Montando a Jacobiana
        for m = 2:6
            %-------------------------------------------------------------------------%
            % Montando os termos de H
            if k == m                                   % Diagonal principal
                H(k-1,m-1) = -V(k)^2 * B(k,k) - V(k)*sum(V.*(G(k,:).*sin(theta(k) - theta) - B(k,:).*cos(theta(k) - theta))); % Montando os termos H da Jacobiana
            else                                        % Resto
                H(k-1,m-1) = V(k)*V(m) * (G(k,m)*sin(theta(k)-theta(m)) - B(k,m)*cos(theta(k) - theta(m))); % Montando os termos H da Jacobiana
            end
            %-------------------------------------------------------------------------%
            % Montando os termos N
            if m >= 4
                if k == m
                    N(k-1,m-3) = V(k) * G(k,k) + sum(V.*(G(k,:).*cos(theta(k) - theta) + B(k,:).*sin(theta(k) - theta)));
                else
                    N(k-1,m-3) = V(k) * (G(k,m) * cos(theta(k) - theta(m)) + B(k,m)*sin(theta(k) - theta(m))); % Montando os termos N da Jacobiana
                end
            end
            %-------------------------------------------------------------------------%
            % Montando os termos M
            if k >= 4
                if k == m
                    M(k-3,m-1) = -V(k)^2 * G(k,k) + V(k)*sum(V.*(G(k,:).*cos(theta(k) - theta) + B(k,:).*sin(theta(k) - theta)));
                else
                    M(k-3,m-1) = -V(k)*V(m) * (G(k,m)*cos(theta(k) - theta(m)) + B(k,m)*sin(theta(k) - theta(m))); % Montando os termos M da Jacobiana
                end
                    
                %-------------------------------------------------------------------------%
                % Montando os termos L
                if m >= 4
                    if k == m
                        L(k-3,m-3) = -V(k) * B(k,k) + sum(V.*(G(k,:).*sin(theta(k) - theta) - B(k,:).*cos(theta(k) - theta))); % Montando os termos L da Jacobiana
                    else
                        L(k-3,m-3) = V(k) * (G(k,m).*sin(theta(k) - theta(m)) - B(k,m).*cos(theta(k) - theta(m))); % Montando os termos L da Jacobiana
                    end
                end
            end
        end       
    end
    J = -[H N; M L];
    
    %calculando novos pontos iniciais
     dx=-inv(J)*gx;  

    %Novos valores de theta
    theta(2:6)=theta(2:6)+dx(1:5)'
    %Novos valores de tensão
    V(4:6)=V(4:6)+dx(6:8)'
  else
        break;
  end;
  
 % Novos valore de gx
    dP2=Pk(2)-V(2)*sum(V.*(G(2,:).*cos(theta(2)-theta) + B(2,:).*sin(theta(2)-theta)));
    dP3=Pk(3)-V(3)*sum(V.*(G(3,:).*cos(theta(3)-theta) + B(3,:).*sin(theta(3)-theta)));
    dP4=Pk(4)-V(4)*sum(V.*(G(4,:).*cos(theta(4)-theta) + B(4,:).*sin(theta(4)-theta)));
    dP5=Pk(5)-V(5)*sum(V.*(G(5,:).*cos(theta(5)-theta) + B(5,:).*sin(theta(5)-theta)));
    dP6=Pk(6)-V(6)*sum(V.*(G(6,:).*cos(theta(6)-theta) + B(6,:).*sin(theta(6)-theta)));
    dQ4=Qk(4)-V(4)*sum(V.*(G(4,:).*sin(theta(4)-theta) - B(4,:).*cos(theta(4)-theta)));
    dQ5=Qk(5)-V(5)*sum(V.*(G(5,:).*sin(theta(5)-theta) - B(5,:).*cos(theta(5)-theta)));
    dQ6=Qk(6)-V(6)*sum(V.*(G(6,:).*sin(theta(6)-theta) - B(6,:).*cos(theta(6)-theta)));
    
    gx=[dP2; dP3; dP4; dP5; dP6; dQ4; dQ5; dQ6];                                    
end


if it==iteracao
    fprintf('-----Sismeta Divergiu-----')
else
    fprintf('Convergiu em %d iterações\n\n',it);
    
    fprintf('Váriaveis Calculadas:\n')
    fprintf('Theta 2 = %f° \n',theta(2)*180/pi());
    fprintf('Theta 3 = %f° \n',theta(3)*180/pi());
    fprintf('Theta 4 = %f° \n',theta(4)*180/pi());
    fprintf('Theta 5 = %f° \n',theta(5)*180/pi());
    fprintf('Theta 6 =  %f° \n',theta(6)*180/pi());
    fprintf('|V4| = %f p.u. \n',V(4));
    fprintf('|V5| = %f p.u. \n',V(5));
    fprintf('|V6| = %f p.u. \n\n',V(6));
end

%% ------- Cálculos para o subsistema 2  -------

if it<iteracao
   Pk(1)= V(1)*sum(V.*(G(1,:).*cos(theta(1)-theta) + B(1,:).*sin(theta(1)-theta)));
   
   Qk(1)=V(1)*sum(V.*(G(1,:).*sin(theta(1)-theta) - B(1,:).*cos(theta(1)-theta)));
   Qk(2)=V(2)*sum(V.*(G(2,:).*sin(theta(2)-theta) - B(2,:).*cos(theta(2)-theta)));
   Qk(3)=V(3)*sum(V.*(G(3,:).*sin(theta(3)-theta) - B(3,:).*cos(theta(3)-theta)));

   
   fprintf('Váriaveis Calculadas:\n')
   fprintf('P1 = %f p.u. \n',Pk(1));
   fprintf('Q1 = %f p.u. \n',Qk(1)); 
   fprintf('Q2 = %f p.u. \n',Qk(2));
   fprintf('Q3 = %f p.u. \n\n',Qk(3));
end

plot(itp,errop) %plotar gráfico do erro
