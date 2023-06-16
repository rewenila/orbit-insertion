function voo_insercao_orbita
% Script para simular um voo de inser��o em �rbita.
% � utilizado o foguete de sondagem VSB30 com os dados abaixo.
close all; clear all; clc;
global Re we mut tg J2 J3 J4 g lc dT Sr fc
global ms m0 mp ti tq ts Isp h0 l_trilho

%% Par�metros de massa estrutural e de carga �til
ms(1) = 284;       % massa estrutural do primeiro est�gio [kg]
ms(2) = 320;       % massa estrutural do segundo est�gio [kg]
% Ajuste dos dois primeiros est�gios
ds(1) = 0.1*ms(1); % redu��o da passa estrutural do primeiro est�gio, convertida em massa de propelente
ds(2) = 0.1*ms(2); % redu��o da passa estrutural do segundo est�gio, convertida em massa de propelente
ms(1) = ms(1) - ds(1);
ms(2) = ms(2) - ds(2);
% Massa estrutural do terceiro est�gio
mo = 400;          % massa original [kg]
ms(3) = 30;        % supor 30kg (valor muito baixo)
% Massa da carga �til
mL = 5;            % supor 5kg (CubeSat 3U mais PPOD)

%% Par�metros propulsivos
Isp(1) = 261;        % impulso espec�fico do primeiro est�gio (motor S31) [s]
Isp(2) = 261;        % impulso espec�fico do segundo est�gio (motor S30) [s]
mp(1) = 677 + ds(1); % massa de propelente do primeiro est�gio (motor S31) [kg], original mais ajuste
mp(2) = 898 + ds(2); % massa de propelente do segundo est�gio (motor S30) [kg], original mais ajuste
ti(1) = 0;           % tempo da igni��o do primeiro est�gio [s]
ti(2) = 15;          % tempo da igni��o do segundo est�gio [s]
tq(1) = 13.5;        % tempo do fim da queima do est�gio 1 (motor S31) [s]
tq(2) = 44;          % tempo do fim da queima do est�gio 2 (motor S30) [s]
ts(1) = 13.5;        % tempo da separa��o do primeiro est�gio [s]
ts(2) = 59;          % tempo da separa��o do segundo est�gio [s]
% Massa de propelente do novo est�gio
mp(3) = mo - ms(3) - mL; % o propelente � tudo que sobra depois da massa de estrutura e de carga �til
% Propelente s�lido com alto impulso espec�fico
Isp(3) = 270; % [s]
% Tempo de igni��o do terceiro est�gio (� um dos par�metros mais
% importantes para se obter a �rbita desejada)
ti(3) = ts(2) + 185; % [s], determinado em rela��o ao tempo de separa��o do segundo est�gio
% Tempo de fim da queima do terceiro est�gio
tq(3) = ti(3) + 30; % [s], determinado por tentativa e erro para gerar a resposta desejada
% Tempo de separa��o do terceiro est�gio
ts(3) = tq(2) + 5; % [s], determinado por chute
% Na mec�nica de voo, este par�metro s� tem import�ncia sobre o arrasto, o
% qual n�o � importante na altitude orbital durante o lan�amento. Na pr�tica,
% est� relacionado a quest�es operacionais tal como a comunica��o e
% monitoramento do sat�lite, visto que o sistema de aquisi��o de dados e
% telemetria do terceiro est�gio � usado para rastrear o sat�lite nos
% primeiros instantes da �rbita.

%% Par�metros aerodin�micos e ambientais
% Fator de corre��o do arrasto a partir de dados de tunel de vento
fc = 1.28;
% �reas de refer�ncia (se��es transversais)
S1 = pi*(0.557/2)^2; % �rea aproximada da se��o transversal do primeiro est�gio [m^2]
S2 = pi*(0.557/2)^2; % �rea aproximada da se��o longitudinal do segundo est�gio [m^2]
S3 = pi*(0.46/2)^2;  % �rea aproximada da se��o longitudinal do terceiro est�gio [m^2]
SL = pi*(0.46/2)^2;  % �rea aproximada da se��o longitudinal da carga �til [m^2]
% Corre��es para levar em conta a �rea molhada, feitas em fun��o do
% comprimento. Assume-se que a corre��o s� ocorre em 50% da �rea de
% refer�ncia, supondo uma influ�ncia de 50% da �rea molhada no arrasto 
% total.
lt = 12.6;       % comprimento total [m]
l2 = lt - 3.214; % comprimento sem o primeiro est�gio [m]
% Comprimento da carga �til
l4 = 0.5; % [m], um chute para um CubeSat 3U, PPOD e acoplamentos
% Comprimento do terceiro est�gio
l3 = l2 - 3.294 - l4; % comprimento do terceiro est�gio [m]
% Fatores de corre��o da �rea molhada (para ajuste da �rea de refer�ncia do
% arrasto)
f2 = (l2/lt)*0.5 + 0.5; % fator de corre��o do segundo est�gio 
f3 = (l3/lt)*0.5 + 0.5; % fator de corre��o do terceiro est�gio
f4 = (l4/lt)*0.5 + 0.5; % fator de corre��o da carga �til
% Vetor de �reas de refer�ncia para c�lculo do arrasto
Sr = [S1
      S2*f2
      S3*f3
      SL*f4];
lc = 0.5; % comprimento caracter�stico
dT = 10;    % Delta-T em rela��o � atmosfera padr�o (15�C no n�vel do mar) [K]

%% Par�metros da Terra - modelo axis sim�trico (WGS-84)
Re = 6378.1370e3;     % raio equatorial da Terra [m]
we = 7.29221150e-5;   % velocidade de rota��o da Terra com respeito ao espa�o inercial [rad/s]
g = 9.80665;          % acelera��o da gravidade padr�o ao n�vel do mar [m/s^2]
mut = 3.986004418e14; % [m3.s^-2]
% Constantes de Jeffery
J2 = 0.00108263;
J3 = -0.00000254;
J4 = -0.00000161;
tg = 0; % tempo em que o meridiano de refer�ncia tem longitude celeste nula [s]

%% Condi��es iniciais - Centro espacial de Alc�ntara
h0 = 0; % altitude da base de lan�amento [m]
delta0 = -2.3267844*pi/180; % latitude inicial [rad]
lon0 = -44.4111042*pi/180;  % longitude inicial [rad]
% Corre��o para compatibilizar com o tra�ado de mapas no MATLAB
lon0 = lon0 - 0.03*pi/180;
% Comprimento do trilho de lan�amento [m]
l_trilho = 10; 

%% Par�metros calculados
% Massa inicial do foguete
m0 = sum(mp) + sum(ms) + mL;
% Dist�ncia radial inicial
r0 = Re + h0;

%% Estudo simplificado pela equa��o de foguete
% Raz�es estruturais
sigma = ms./(ms + mp);
% Massa total na decolagem
m01 = m0;
% Massa total na igni��o do segundo est�gio
m02 = m01 - mp(1) - ms(1);
% Massa inicial antes da queima do terceiro est�gio
m03 = m02 - mp(2) - ms(2);
% Raz�o de carga �til do primeiro e segundo est�gios
lamb(1) = m02/m01;
lamb(2) = m03/m02;
% Raz�o de carga �til do terceiro est�gio
lamb(3) = mL/m03;
% Raz�o de carga �til total
lambL = prod(lamb);
% Velocidade de exaust�o
ve = g*Isp;
% Delta v total ideal
Dv = -sum(ve.*log(sigma + (1 - sigma).*lamb));
% Delta v ideal impresso pelo terceiro est�gio
Dv3 = -ve(3)*log(sigma(3) + (1 - sigma(3))*lamb(3));

%% Mostra dados na tela
disp('�rea de refer�ncia do foguete com primeiro est�gio (m^2):');
disp(Sr(1));
disp('�rea de refer�ncia do foguete com segundo est�gio (m^2):');
disp(Sr(2));
disp('�rea de refer�ncia do foguete com terceiro est�gio (m^2):');
disp(Sr(3));
disp('�rea de refer�ncia da arga �til (m^2):');
disp(Sr(4));
disp('Massa inicial antes da queima de cada est�gio (kg):');
disp(m0);
disp('Raz�es estruturais:');
disp(sigma);
disp('Raz�es de carga �til:');
disp(lamb);
disp('Raz�o de carga �til total:');
disp(lambL);
disp('Velocidades de exaust�o (m/s):');
disp(ve);
disp('Impulso de velocidade total ideal (m/s):');
disp(Dv);
disp('Impulso de velocidade ideal impresso pelo terceiro est�gio (m/s):');
disp(Dv3);

%% Verifica se o usu�rio deseja continual
cont = input('Deseja continuar? (1 = sim) (0 = n�o): ');
if cont == 0
    return
end

%% Simula��o
simula = 1;
while simula == 1
    %% Dados fornecidos pelo usu�rio
    TF = input('Informe o tempo de simula��o (s): ');
    v0 = input('Informe o valor inicial da velocidade relativa (m/s): ');
    phi0 = input('Informe a condi��o inicial do �ngulo de eleva��o (graus): ');
    phi0 = phi0*pi/180;
    in = input('Informe a inclina��o da �rbita desejada (graus): ');
    in = in*pi/180;
    % Faz um teste de factibilidade usando a inclina��o da �rbita e a
    % latitude inicial
    y = cos(in)/cos(delta0);
    if abs(y) > 1
        disp('N�o � poss�vel atingir a inclina��o a partir da latitude inicial. Calculando a maior poss�vel.');
        y = sign(y);
    end
    Ai_f = asin(y); % condi��o final do azimute de velocidade inercial
    % Estimativa da condi��o inicial de azimute de velocidade relativa
    % Raio de uma �rbita circular de refer�ncia com 100km de altitude
    rref = Re + 100e3;
    % Velocidade de uma �rbita de refer�ncia com 100km de altitude
    viref = sqrt(mut/rref);
    A0 = atan(tan(Ai_f) - (rref*we*cos(delta0))/(viref*cos(Ai_f)));
    disp('Condi��o final de azimute de velocidade inercial (�): ');
    disp(Ai_f*180/pi);
    disp('Condi��o inicial de azimute de velocidade relativa (�): ');
    disp(A0*180/ pi);
    
    % Condi��o inicial
    X0 = [v0;A0;phi0;r0;delta0;lon0];
    %options = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',0.5);
    options = odeset('RelTol',1e-8,'AbsTol',1e-12);
    %[t,X] = ode45('dinamica_foguete',[0 TF],X0,options);
    [t,X] = ode15s('dinamica_foguete',[0 TF],X0,options);

    %% P�s processamento
    N = length(t);      % n�mero de instantes de tempo
    V = zeros(N,1); A = zeros(N,1); phi = zeros(N,1); % magnitude, azimute e eleva��o da velocidade relativa
    h = zeros(N,1); delta = zeros(N,1); lon = zeros(N,1); % altitude, latitude e longitude do referencial fixo ao planeta
    m = zeros(N,1);     % massa
    ft = zeros(N,1);    % for�a propulsiva
    D = zeros(N,1);     % for�a de arrasto
    q = zeros(N,1);     % press�o din�mica
    M = zeros(N,1);     % n�mero de Mach
    T = zeros(N,1);     % temperatura
    rho = zeros(N,1);   % densidade
    Vi = zeros(N,1);    % magnitude da velocida inercial
    phii = zeros(N,1);  % eleva��o da velocidade inercial
    Ai = zeros(N,1);    % azimute da velocidade inercial
    longc = zeros(N,1); % longitude celeste
    ee = zeros(N,1);    % energia espec�fica
    a = zeros(N,1);     % semi eixo maior da �rbita
    e = zeros(N,1);     % excentricidade da �rbita
    tau = zeros(N,1);   % tempo de periastro
    OM = zeros(N,1);    % ascen��o reta do nodo ascendente
    in = zeros(N,1);    % inclina��o da �rbita
    om = zeros(N,1);    % argumento de periastro
    R0 = zeros(N,3);    % posi��o no referencial ECI 

    for i = 1:N
        V(i) = X(i,1);      % velocidade relativa
        A(i) = X(i,2);      % azimute da velocidade reltiva
        phi(i) = X(i,3);    % eleva��o da velocidade relativa
        h(i) = X(i,4) - Re; % altitude
        r = X(i,4);         % dist�ncia radial
        delta(i) = X(i,5);  % latitude
        lon(i) = X(i,6);    % longitude no referencial fixo ao planeta
        [ft(i),m(i)] = propulsao_N_estagios(t(i)); % for�a propulsiva e massa
        [T(i),~,rho(i),~,M(i),~,~,Kn,~,~] = atm_padrao(h(i),V(i),lc,dT); % par�metros atmosf�ricos
        [D(i),~,~] = aerodinamica_N_estagios(t(i),V(i),h(i),M(i),Kn,T(i),rho(i)); % for�as aerodin�micas
        q(i) = 0.5*rho(i)*V(i)^2; % press�o din�mica
        % Coordenadas da velocidade inercial no referencial horizontal local
        [Vi(i),phii(i),Ai(i)] = Vrel2Vine(V(i),phi(i),A(i),we,r,delta(i));
        % Longitude celeste 
        longc(i) = long_ECEF2ECI(t(i),lon(i),we,tg);
        % Energia espec�fica da �rbita
        ee(i) = Vi(i)^2/2 - mut/r;
        % Posi��o e velocidade inercial no referencial celeste ECI
        [rc0,vc0] = RvelPolar2RvelRet(Vi(i),Ai(i),phii(i),r,delta(i),longc(i));
        R0(i,:) = rc0';
        % Par�metros orbitais
        par_orb = det_orb(t(i),rc0,vc0,mut);
        a(i) = par_orb(1); e(i) = par_orb(2); tau(i) = par_orb(3); 
        OM(i) = par_orb(4); in(i) = par_orb(5); om(i) = par_orb(6); 
    end
    
    %% An�lise de �rbita
    % Gera��o de vetores para tra�ar gr�fico
    % Semi eixo maior de uma �rbita circular de refer�ncia com 100 km de
    % altitude 
    ar = rref*ones(1,N);
    % Velocidade de uma �rbita de refer�ncia com 100 km de altitude
    Vir = viref*ones(1,N);
    % Energia espec�fica de uma �rbita de refer�ncia com 100 km de altitude
    eer = -mut./(2*ar);
    % Altitude e velocidade inercial no fim da queima do terceiro est�gio
    for i = 1:N
        if t(i) > tq(3)
            break;
        end
    end
    ifq = i - 1;
    tfq = t(ifq); % tempo do fim da queima do terceiro est�gio
    Vfq = Vi(ifq)*ones(1,N); % velocidade inercial no fim da queima do terceiro est�gio
    hfq = h(ifq)*ones(1,N); % altitude no fim da queima do terceiro est�gio
    P = 2*pi*sqrt((Re + hfq(1))^3/mut); % per�odo da �rbita obtida
    disp('Per�odo da �rbita obtida (min): ');
    disp(P/60);
    rp = a(ifq)*(1 - e(ifq)); % raio do periastro
    ra = a(ifq)*(1 + e(ifq)); % raio do apoastro
    disp('Raio do periastro (km): ');
    disp(rp/1e3);
    disp('Raio do apoastro (km): ');
    disp(ra/1e3);
    disp('Altitude do periastro (km): ');
    disp((rp - Re)/1e3);
    disp('Altitude do apoastro (km): ');
    disp((ra - Re)/1e3);
    
    %% Gr�ficos
    figure(1)
    subplot(231)
    plot(t,V,'LineWidth',2); grid; axis tight; xlabel('t(s)'); ylabel('V(m/s)');
    subplot(232)
    plot(t,A*180/pi,'LineWidth',2); grid; axis tight; xlabel('t(s)'); ylabel('A(�)');
    subplot(233)
    plot(t,phi*180/pi,'LineWidth',2); grid; axis tight; xlabel('t(s)'); ylabel('\phi(�)');
    subplot(234)
    plot(t,h/1e3,'LineWidth',2); hold;
    plot(t,hfq/1e3,'--',tfq,hfq(1)/1e3,'*');
    grid; axis tight; xlabel('t(s)'); ylabel('h(km)');
    legend('altitude','altitude no fim da queima do 3� est�gio');
    subplot(235)
    plot(t,delta*180/pi,'LineWidth',2); grid; axis tight; xlabel('t(s)'); ylabel('\delta(�)');
    subplot(236)
    plot(t,lon*180/pi,'LineWidth',2); grid; axis tight; xlabel('t(s)'); ylabel('l(�)');

    figure(2)
    subplot(221); 
    plot(t,Vi,'LineWidth',2); hold;
    plot(t,Vir,'--',t,Vfq,'-.',tfq,Vfq(1),'*'); grid; xlabel('t(s)'); ylabel('V_i(m/s)');
    legend('Velocidade inercial','Velocidade de �rbita circular de refer�ncia de 100km de altitude',...
        'Velocidade no fim da queima do terceiro est�gio');
    subplot(222);
    plot(t,Ai*180/pi,'LineWidth',2); grid; axis tight; xlabel('t(s)'); ylabel('A_i(�)');
    subplot(223);
    plot(t,phii*180/pi,'LineWidth',2); grid; axis tight; xlabel('t(s)'); ylabel('\phi_i(�)');
    subplot(224);
    plot(t,longc*180/pi,'LineWidth',2); grid; axis tight; xlabel('t(s)'); ylabel('\lambda(�)');   
    
    figure(3)
    subplot(211)
    plot(t,ft,'LineWidth',2); grid; axis tight; xlabel('t(s)'); ylabel('f_t(N)');
    subplot(212)
    plot(t,m,'LineWidth',2); grid; axis tight; xlabel('t(s)'); ylabel('m(kg)');

    figure(4)
    subplot(311)
    plot(t,D,'LineWidth',2); grid; axis tight; xlabel('t(s)'); ylabel('D(N)');
    subplot(323)
    plot(t,q,'LineWidth',2); grid; axis tight; xlabel('t(s)'); ylabel('q(N/m^2)');
    subplot(324)
    plot(t,M,'LineWidth',2); grid; axis tight; xlabel('t(s)'); ylabel('M(-)');
    subplot(325)
    plot(t,T-273.15,'LineWidth',2); grid; axis tight; xlabel('t(s)'); ylabel('T(�C)');
    subplot(326)
    plot(t,rho,'LineWidth',2); grid; axis tight; xlabel('t(s)'); ylabel('\rho(kg/m^3)');

    figure(5)
    subplot(331)
    plot(t,a/1e3,'LineWidth',2); hold; 
    legend('Semi eixo maior');
    plot(t,ar/1e3,'--',t,Re*ones(1,N)/1e3,'-.');
    grid; xlabel('t(s)'); ylabel('a(km)');
    legend('Semi eixo maior de �rbita com 100km','Raio da Terra');
    subplot(332)
    plot(t,e,'LineWidth',2); grid; axis tight; xlabel('t(s)'); ylabel('e(-)');
    subplot(333)
    plot(t,tau,'LineWidth',2); grid; axis tight; xlabel('t(s)'); ylabel('\tau(s)');
    subplot(334)
    plot(t,OM*180/pi,'LineWidth',2); grid; axis tight; xlabel('t(s)'); ylabel('\Omega(�)');
    subplot(335)
    plot(t,in*180/pi,'LineWidth',2); grid; axis tight; xlabel('t(s)'); ylabel('i(�)');
    subplot(336)
    plot(t,om*180/pi,'LineWidth',2); grid; axis tight; xlabel('t(s)'); ylabel('\omega(�)');
    subplot(325)
    plot(t,ee,t,eer,'--','LineWidth',2); grid; axis tight; xlabel('t(s)'); ylabel('\epsilon(J/kg)');
    legend('Energia espec�fica','Energia espec�fica de �rbita com 100km');
    
    figure(6)
    traj = [delta lon]*180/pi;
    desenha_mapa_trajetoria([delta0*180/pi,lon0*180/pi,h0],traj);
    
    figure(7)
    plot3(R0(:,1)/1e3,R0(:,2)/1e3,R0(:,3)/1e3,'LineWidth',2);
    xlabel('X(km)'); ylabel('Y(km)'); zlabel('Z(km)'); axis tight;
    hold on
    h1 = gca;
    earth_sphere(h1,'km');
    
    simula = input('Deseja simular novamente? (1 = sim) (0 = n�o): ');

end
    
end











