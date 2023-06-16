function [T,p,rho,ainf,M,mu,Pr,Kn,d,Re] = atm_padrao(h,v,lc,dT)
%% Fun��o para c�lculo da atmosfera padr�o para a altitude geom�trica de 0 a 2000km
% Entradas
% h: altitude geom�trica [m]
% v: velocidade do ve�culo em rela��o ao escoamento n�o perturbado [m/s]
% lc: comprimento caracter�stico do ve�culo [m]
% Sa�das
% T: temperatura [K]
% p: press�o [Pa = N/m^2]
% rho: densidade [kg/m^3]
% ainf: velocidade do som [m/s]
% M: n�mero de Mach [adm]
% mu: coeficiente de viscosidade din�mica [kg/m*s]
% Pr: n�mero de Prandt [adm]
% Kn: n�mero de Knudsen [adm]
% d: par�metro de regime de escoamento [adm]
% Re: n�mero de Reynolds [adm]

%% Entrada de dados
% Vetores com os seguintes dados, sendo que cada linha � uma camada do
% modelo de atmosfera padr�o: altitude no in�cio da camada [m], temperatura
% no in�cio da camada [K], constante de g�s ideal do ar na camada [J/kg.K],
% taxa de lapso t�rmico [K/m]
hi = 1e3*[0
          11.0191
          20.0631
          32.1619
          47.3501
          51.4125
          71.8020
          86
          100
          110
          120
          150
          160
          170
          190
          230
          300
          400
          500
          600
          700];
      
Ti = [288.15
      216.65
      216.65
      228.65
      270.65
      270.65
      214.65
      186.946
      210.02
      257.0
      349.49
      892.79
      1022.2
      1103.4
      1205.4
      1322.3
      1423.1
      1487.4
      1506.1
      1506.1
      1507.6];
  
R = [287.0
     287.0
     287.0
     287.0
     287.0
     287.0
     287.0
     287.02
     287.02
     287.84
     291.06
     308.79
     311.80
     313.69
     321.57
     336.68
     366.84
     416.88
     463.36
     493.63
     514.08
     514.08];
 
 a = 1e-3*[-6.5
           0.0
           1.0
           2.8
           0.0
           -2.8
           -2.0
           1.693
           5.0
           10.0
           20.0
           15.0
           10.0
           7.0
           5.0
           4.0
           3.3
           2.6
           1.7
           1.1
           0.0];
 
% Constantes
g0 = 9.80665; % valor ao n�vel do mar da acelera��o da gravidade [m/s^2]
Na = 6.0220978e23; % n�mero de Avogadro
sigma = 3.65e-10; % di�metro de colis�o para o ar [m]
m0 = 28.964e-3; % massa molar ao n�vel do mar [kg/Mol]
P0 = 1.01325e5; % press�o padr�o ao n�vel do mar [N/m^2]
Re = 6378.14e3; % raio m�dio da Terra [m]
gamma = 1.405; % raz�o de calores espec�ficos ao n�vel do mar

%% Constantes calculadas
% N�mero beta associado � dist�ncia radial m�dia ao n�vel do mar
beta = 2/Re;

%% Identifica a camada � qual a altitude pertence
% Contador
i = 1;
if h < 0
    disp('Cuidado: altitude negativa.');
    disp('Os resultados apresentados dizem respeito a h = 0.');
    i = 1; h = 0;
elseif h > 2000e3
    disp('A altitude fornecida est� acima do limite superior de 2000 km.');
    disp('Os resultados dizem respeito a h = 2000 km.');
    i = 21; h = 2000e3;
else
    for i = 1:21
        if (i == 21)
            break
        elseif ((h >= hi(i)) && (h < hi(i+1)))
            break
        end
    end
end

%% Realiza os c�lculos
% Press�o no in�cio da camada i
Pi = P0; % press�o inicial da camada, inicializa com o valor ao n�vel do mar
px = P0; % vari�vel auxiliar para guardar o valor de press�o inicial na camada anterior
for j = 2:i
    if a(j-1) ~= 0 % verifica se a camada n�o � isot�rmica
        % Calcula a press�o inicial da camada i a partir do modelo da
        % camada i-1
        A = 1 + (a(j-1)*(hi(j) - hi(j-1)))/Ti(j-1);
        B = -(g0/(R(j-1)*a(j-1)))*(1 + beta*(Ti(j-1)/a(j-1) - hi(j-1)));
        C = (g0*beta/(R(j-1)*a(j-1)))*(hi(j) - hi(j-1));
        Pi = px*(A^B)*exp(C);
        px = Pi; % o valor atual ser� o valor anterior na pr�xima itera��o
    else 
        % Calcula a press�o inicial da camada i pelo modelo isot�rmico
        Pi = px*exp(-(g0/(R(j-1)*Ti(j-1)))*(hi(j) - hi(j-1))*(1 - (beta/2)*...
            (hi(j) - hi(j-1))));
        px = Pi; % o valor atual ser� o valor anterior na pr�xima itera��o
    end
end

% Temperatura padr�o
T = Ti(i) + a(i)*(h - hi(i));
% Corrige pelo valor de delta T adotado
T = T + dT;
% Press�o
if a(i) ~= 0 % verifica se a camada n�o � isot�rmica
    % Calcula a press�o
    A = 1 + (a(i)*(h - hi(i)))/Ti(i);
    B = -(g0/(R(i)*a(i)))*(1 + beta*(Ti(i)/a(i) - hi(i)));
    C = (g0*beta/(R(i)*a(i)))*(h - hi(i));
    p = Pi*(A^B)*exp(C);
else
    % Calcula a press�o pelo modelo isot�rmico
    p = Pi*exp(-(g0/(R(i)*Ti(i)))*(h - hi(i))*(1 - (beta/2)*(h - hi(i))));
end
% Densidade
rho = p/(R(i)*T);
% Velocidade do som
ainf = sqrt(gamma*R(i)*T);
% N�mero de Mach
M = v/ainf;
% Coeficiente de viscosidade din�mica
mu = 1.458e-6*(T^(3/2))/(T + 110.4);
% N�mero de Prandt
cp = R(i)*gamma/(gamma - 1);
kT = (2.64638e-3*(T^(3/2)))/(T + 245.4*(10^(-12/T)));
Pr = mu*cp/kT;
% N�mero de Knudsen
lam = m0/(sqrt(2)*pi*sigma^2*rho*Na);
Kn = lam/lc;
% Par�metro de regime de escoamento
if Kn >= 10
    d = 1;
elseif Kn <= 0.01
    d = 2;
else
    d = 3;
end
% N�mero de Reynolds
Re = rho*v*lc/mu;

end
