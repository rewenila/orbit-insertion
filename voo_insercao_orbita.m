function voo_insercao_orbita
% Script para simular um voo de inserção em órbita.
% É utilizado o foguete de sondagem VSB30 com os dados abaixo.
close all; clear all; clc;
global Re we mut tg J2 J3 J4 g lc dT Sr fc
global ms m0 mp ti tq ts Isp h0 l_trilho

%% Parâmetros de massa estrutural e de carga útil
ms(1) = 284;       % massa estrutural do primeiro estágio [kg]
ms(2) = 320;       % massa estrutural do segundo estágio [kg]
% Ajuste dos dois primeiros estágios
ds(1) = 0.1*ms(1); % redução da passa estrutural do primeiro estágio, convertida em massa de propelente
ds(2) = 0.1*ms(2); % redução da passa estrutural do segundo estágio, convertida em massa de propelente
ms(1) = ms(1) - ds(1);
ms(2) = ms(2) - ds(2);
% Massa estrutural do terceiro estágio
mo = 400;          % massa original [kg]
ms(3) = 30;        % supor 30kg (valor muito baixo)
% Massa da carga útil
mL = 5;            % supor 5kg (CubeSat 3U mais PPOD)

%% Parâmetros propulsivos
Isp(1) = 261;        % impulso específico do primeiro estágio (motor S31) [s]
Isp(2) = 261;        % impulso específico do segundo estágio (motor S30) [s]
mp(1) = 677 + ds(1); % massa de propelente do primeiro estágio (motor S31) [kg], original mais ajuste
mp(2) = 898 + ds(2); % massa de propelente do segundo estágio (motor S30) [kg], original mais ajuste
ti(1) = 0;           % tempo da ignição do primeiro estágio [s]
ti(2) = 15;          % tempo da ignição do segundo estágio [s]
tq(1) = 13.5;        % tempo do fim da queima do estágio 1 (motor S31) [s]
tq(2) = 44;          % tempo do fim da queima do estágio 2 (motor S30) [s]
ts(1) = 13.5;        % tempo da separação do primeiro estágio [s]
ts(2) = 59;          % tempo da separação do segundo estágio [s]
% Massa de propelente do novo estágio
mp(3) = mo - ms(3) - mL; % o propelente é tudo que sobra depois da massa de estrutura e de carga útil
% Propelente sólido com alto impulso específico
Isp(3) = 270; % [s]
% Tempo de ignição do terceiro estágio (é um dos parâmetros mais
% importantes para se obter a órbita desejada)
ti(3) = ts(2) + 185; % [s], determinado em relação ao tempo de separação do segundo estágio
% Tempo de fim da queima do terceiro estágio
tq(3) = ti(3) + 30; % [s], determinado por tentativa e erro para gerar a resposta desejada
% Tempo de separação do terceiro estágio
ts(3) = tq(2) + 5; % [s], determinado por chute
% Na mecânica de voo, este parâmetro só tem importância sobre o arrasto, o
% qual não é importante na altitude orbital durante o lançamento. Na prática,
% está relacionado a questões operacionais tal como a comunicação e
% monitoramento do satélite, visto que o sistema de aquisição de dados e
% telemetria do terceiro estágio é usado para rastrear o satélite nos
% primeiros instantes da órbita.

%% Parâmetros aerodinâmicos e ambientais
% Fator de correção do arrasto a partir de dados de tunel de vento
fc = 1.28;
% Áreas de referência (seções transversais)
S1 = pi*(0.557/2)^2; % área aproximada da seção transversal do primeiro estágio [m^2]
S2 = pi*(0.557/2)^2; % área aproximada da seção longitudinal do segundo estágio [m^2]
S3 = pi*(0.46/2)^2;  % área aproximada da seção longitudinal do terceiro estágio [m^2]
SL = pi*(0.46/2)^2;  % área aproximada da seção longitudinal da carga útil [m^2]
% Correções para levar em conta a área molhada, feitas em função do
% comprimento. Assume-se que a correção só ocorre em 50% da área de
% referência, supondo uma influência de 50% da área molhada no arrasto 
% total.
lt = 12.6;       % comprimento total [m]
l2 = lt - 3.214; % comprimento sem o primeiro estágio [m]
% Comprimento da carga útil
l4 = 0.5; % [m], um chute para um CubeSat 3U, PPOD e acoplamentos
% Comprimento do terceiro estágio
l3 = l2 - 3.294 - l4; % comprimento do terceiro estágio [m]
% Fatores de correção da área molhada (para ajuste da área de referência do
% arrasto)
f2 = (l2/lt)*0.5 + 0.5; % fator de correção do segundo estágio 
f3 = (l3/lt)*0.5 + 0.5; % fator de correção do terceiro estágio
f4 = (l4/lt)*0.5 + 0.5; % fator de correção da carga útil
% Vetor de áreas de referência para cálculo do arrasto
Sr = [S1
      S2*f2
      S3*f3
      SL*f4];
lc = 0.5; % comprimento característico
dT = 10;    % Delta-T em relação à atmosfera padrão (15°C no nível do mar) [K]

%% Parâmetros da Terra - modelo axis simétrico (WGS-84)
Re = 6378.1370e3;     % raio equatorial da Terra [m]
we = 7.29221150e-5;   % velocidade de rotação da Terra com respeito ao espaço inercial [rad/s]
g = 9.80665;          % aceleração da gravidade padrão ao nível do mar [m/s^2]
mut = 3.986004418e14; % [m3.s^-2]
% Constantes de Jeffery
J2 = 0.00108263;
J3 = -0.00000254;
J4 = -0.00000161;
tg = 0; % tempo em que o meridiano de referência tem longitude celeste nula [s]

%% Condições iniciais - Centro espacial de Alcântara
h0 = 0; % altitude da base de lançamento [m]
delta0 = -2.3267844*pi/180; % latitude inicial [rad]
lon0 = -44.4111042*pi/180;  % longitude inicial [rad]
% Correção para compatibilizar com o traçado de mapas no MATLAB
lon0 = lon0 - 0.03*pi/180;
% Comprimento do trilho de lançamento [m]
l_trilho = 10; 

%% Parâmetros calculados
% Massa inicial do foguete
m0 = sum(mp) + sum(ms) + mL;
% Distância radial inicial
r0 = Re + h0;

%% Estudo simplificado pela equação de foguete
% Razões estruturais
sigma = ms./(ms + mp);
% Massa total na decolagem
m01 = m0;
% Massa total na ignição do segundo estágio
m02 = m01 - mp(1) - ms(1);
% Massa inicial antes da queima do terceiro estágio
m03 = m02 - mp(2) - ms(2);
% Razão de carga útil do primeiro e segundo estágios
lamb(1) = m02/m01;
lamb(2) = m03/m02;
% Razão de carga útil do terceiro estágio
lamb(3) = mL/m03;
% Razão de carga útil total
lambL = prod(lamb);
% Velocidade de exaustão
ve = g*Isp;
% Delta v total ideal
Dv = -sum(ve.*log(sigma + (1 - sigma).*lamb));
% Delta v ideal impresso pelo terceiro estágio
Dv3 = -ve(3)*log(sigma(3) + (1 - sigma(3))*lamb(3));

%% Mostra dados na tela
disp('Área de referência do foguete com primeiro estágio (m^2):');
disp(Sr(1));
disp('Área de referência do foguete com segundo estágio (m^2):');
disp(Sr(2));
disp('Área de referência do foguete com terceiro estágio (m^2):');
disp(Sr(3));
disp('Área de referência da arga útil (m^2):');
disp(Sr(4));
disp('Massa inicial antes da queima de cada estágio (kg):');
disp(m0);
disp('Razões estruturais:');
disp(sigma);
disp('Razões de carga útil:');
disp(lamb);
disp('Razão de carga útil total:');
disp(lambL);
disp('Velocidades de exaustão (m/s):');
disp(ve);
disp('Impulso de velocidade total ideal (m/s):');
disp(Dv);
disp('Impulso de velocidade ideal impresso pelo terceiro estágio (m/s):');
disp(Dv3);

%% Verifica se o usuário deseja continual
cont = input('Deseja continuar? (1 = sim) (0 = não): ');
if cont == 0
    return
end

%% Simulação
simula = 1;
while simula == 1
    %% Dados fornecidos pelo usuário
    TF = input('Informe o tempo de simulação (s): ');
    v0 = input('Informe o valor inicial da velocidade relativa (m/s): ');
    phi0 = input('Informe a condição inicial do ângulo de elevação (graus): ');
    phi0 = phi0*pi/180;
    in = input('Informe a inclinação da órbita desejada (graus): ');
    in = in*pi/180;
    % Faz um teste de factibilidade usando a inclinação da órbita e a
    % latitude inicial
    y = cos(in)/cos(delta0);
    if abs(y) > 1
        disp('Não é possível atingir a inclinação a partir da latitude inicial. Calculando a maior possível.');
        y = sign(y);
    end
    Ai_f = asin(y); % condição final do azimute de velocidade inercial
    % Estimativa da condição inicial de azimute de velocidade relativa
    % Raio de uma órbita circular de referência com 100km de altitude
    rref = Re + 100e3;
    % Velocidade de uma órbita de referência com 100km de altitude
    viref = sqrt(mut/rref);
    A0 = atan(tan(Ai_f) - (rref*we*cos(delta0))/(viref*cos(Ai_f)));
    disp('Condição final de azimute de velocidade inercial (°): ');
    disp(Ai_f*180/pi);
    disp('Condição inicial de azimute de velocidade relativa (°): ');
    disp(A0*180/ pi);
    
    % Condição inicial
    X0 = [v0;A0;phi0;r0;delta0;lon0];
    %options = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',0.5);
    options = odeset('RelTol',1e-8,'AbsTol',1e-12);
    %[t,X] = ode45('dinamica_foguete',[0 TF],X0,options);
    [t,X] = ode15s('dinamica_foguete',[0 TF],X0,options);

    %% Pós processamento
    N = length(t);      % número de instantes de tempo
    V = zeros(N,1); A = zeros(N,1); phi = zeros(N,1); % magnitude, azimute e elevação da velocidade relativa
    h = zeros(N,1); delta = zeros(N,1); lon = zeros(N,1); % altitude, latitude e longitude do referencial fixo ao planeta
    m = zeros(N,1);     % massa
    ft = zeros(N,1);    % força propulsiva
    D = zeros(N,1);     % força de arrasto
    q = zeros(N,1);     % pressão dinâmica
    M = zeros(N,1);     % número de Mach
    T = zeros(N,1);     % temperatura
    rho = zeros(N,1);   % densidade
    Vi = zeros(N,1);    % magnitude da velocida inercial
    phii = zeros(N,1);  % elevação da velocidade inercial
    Ai = zeros(N,1);    % azimute da velocidade inercial
    longc = zeros(N,1); % longitude celeste
    ee = zeros(N,1);    % energia específica
    a = zeros(N,1);     % semi eixo maior da órbita
    e = zeros(N,1);     % excentricidade da órbita
    tau = zeros(N,1);   % tempo de periastro
    OM = zeros(N,1);    % ascenção reta do nodo ascendente
    in = zeros(N,1);    % inclinação da órbita
    om = zeros(N,1);    % argumento de periastro
    R0 = zeros(N,3);    % posição no referencial ECI 

    for i = 1:N
        V(i) = X(i,1);      % velocidade relativa
        A(i) = X(i,2);      % azimute da velocidade reltiva
        phi(i) = X(i,3);    % elevação da velocidade relativa
        h(i) = X(i,4) - Re; % altitude
        r = X(i,4);         % distância radial
        delta(i) = X(i,5);  % latitude
        lon(i) = X(i,6);    % longitude no referencial fixo ao planeta
        [ft(i),m(i)] = propulsao_N_estagios(t(i)); % força propulsiva e massa
        [T(i),~,rho(i),~,M(i),~,~,Kn,~,~] = atm_padrao(h(i),V(i),lc,dT); % parâmetros atmosféricos
        [D(i),~,~] = aerodinamica_N_estagios(t(i),V(i),h(i),M(i),Kn,T(i),rho(i)); % forças aerodinâmicas
        q(i) = 0.5*rho(i)*V(i)^2; % pressão dinâmica
        % Coordenadas da velocidade inercial no referencial horizontal local
        [Vi(i),phii(i),Ai(i)] = Vrel2Vine(V(i),phi(i),A(i),we,r,delta(i));
        % Longitude celeste 
        longc(i) = long_ECEF2ECI(t(i),lon(i),we,tg);
        % Energia específica da órbita
        ee(i) = Vi(i)^2/2 - mut/r;
        % Posição e velocidade inercial no referencial celeste ECI
        [rc0,vc0] = RvelPolar2RvelRet(Vi(i),Ai(i),phii(i),r,delta(i),longc(i));
        R0(i,:) = rc0';
        % Parâmetros orbitais
        par_orb = det_orb(t(i),rc0,vc0,mut);
        a(i) = par_orb(1); e(i) = par_orb(2); tau(i) = par_orb(3); 
        OM(i) = par_orb(4); in(i) = par_orb(5); om(i) = par_orb(6); 
    end
    
    %% Análise de órbita
    % Geração de vetores para traçar gráfico
    % Semi eixo maior de uma órbita circular de referência com 100 km de
    % altitude 
    ar = rref*ones(1,N);
    % Velocidade de uma órbita de referência com 100 km de altitude
    Vir = viref*ones(1,N);
    % Energia específica de uma órbita de referência com 100 km de altitude
    eer = -mut./(2*ar);
    % Altitude e velocidade inercial no fim da queima do terceiro estágio
    for i = 1:N
        if t(i) > tq(3)
            break;
        end
    end
    ifq = i - 1;
    tfq = t(ifq); % tempo do fim da queima do terceiro estágio
    Vfq = Vi(ifq)*ones(1,N); % velocidade inercial no fim da queima do terceiro estágio
    hfq = h(ifq)*ones(1,N); % altitude no fim da queima do terceiro estágio
    P = 2*pi*sqrt((Re + hfq(1))^3/mut); % período da órbita obtida
    disp('Período da órbita obtida (min): ');
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
    
    %% Gráficos
    figure(1)
    subplot(231)
    plot(t,V,'LineWidth',2); grid; axis tight; xlabel('t(s)'); ylabel('V(m/s)');
    subplot(232)
    plot(t,A*180/pi,'LineWidth',2); grid; axis tight; xlabel('t(s)'); ylabel('A(°)');
    subplot(233)
    plot(t,phi*180/pi,'LineWidth',2); grid; axis tight; xlabel('t(s)'); ylabel('\phi(°)');
    subplot(234)
    plot(t,h/1e3,'LineWidth',2); hold;
    plot(t,hfq/1e3,'--',tfq,hfq(1)/1e3,'*');
    grid; axis tight; xlabel('t(s)'); ylabel('h(km)');
    legend('altitude','altitude no fim da queima do 3° estágio');
    subplot(235)
    plot(t,delta*180/pi,'LineWidth',2); grid; axis tight; xlabel('t(s)'); ylabel('\delta(°)');
    subplot(236)
    plot(t,lon*180/pi,'LineWidth',2); grid; axis tight; xlabel('t(s)'); ylabel('l(°)');

    figure(2)
    subplot(221); 
    plot(t,Vi,'LineWidth',2); hold;
    plot(t,Vir,'--',t,Vfq,'-.',tfq,Vfq(1),'*'); grid; xlabel('t(s)'); ylabel('V_i(m/s)');
    legend('Velocidade inercial','Velocidade de órbita circular de referência de 100km de altitude',...
        'Velocidade no fim da queima do terceiro estágio');
    subplot(222);
    plot(t,Ai*180/pi,'LineWidth',2); grid; axis tight; xlabel('t(s)'); ylabel('A_i(°)');
    subplot(223);
    plot(t,phii*180/pi,'LineWidth',2); grid; axis tight; xlabel('t(s)'); ylabel('\phi_i(°)');
    subplot(224);
    plot(t,longc*180/pi,'LineWidth',2); grid; axis tight; xlabel('t(s)'); ylabel('\lambda(°)');   
    
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
    plot(t,T-273.15,'LineWidth',2); grid; axis tight; xlabel('t(s)'); ylabel('T(°C)');
    subplot(326)
    plot(t,rho,'LineWidth',2); grid; axis tight; xlabel('t(s)'); ylabel('\rho(kg/m^3)');

    figure(5)
    subplot(331)
    plot(t,a/1e3,'LineWidth',2); hold; 
    legend('Semi eixo maior');
    plot(t,ar/1e3,'--',t,Re*ones(1,N)/1e3,'-.');
    grid; xlabel('t(s)'); ylabel('a(km)');
    legend('Semi eixo maior de órbita com 100km','Raio da Terra');
    subplot(332)
    plot(t,e,'LineWidth',2); grid; axis tight; xlabel('t(s)'); ylabel('e(-)');
    subplot(333)
    plot(t,tau,'LineWidth',2); grid; axis tight; xlabel('t(s)'); ylabel('\tau(s)');
    subplot(334)
    plot(t,OM*180/pi,'LineWidth',2); grid; axis tight; xlabel('t(s)'); ylabel('\Omega(°)');
    subplot(335)
    plot(t,in*180/pi,'LineWidth',2); grid; axis tight; xlabel('t(s)'); ylabel('i(°)');
    subplot(336)
    plot(t,om*180/pi,'LineWidth',2); grid; axis tight; xlabel('t(s)'); ylabel('\omega(°)');
    subplot(325)
    plot(t,ee,t,eer,'--','LineWidth',2); grid; axis tight; xlabel('t(s)'); ylabel('\epsilon(J/kg)');
    legend('Energia específica','Energia específica de órbita com 100km');
    
    figure(6)
    traj = [delta lon]*180/pi;
    desenha_mapa_trajetoria([delta0*180/pi,lon0*180/pi,h0],traj);
    
    figure(7)
    plot3(R0(:,1)/1e3,R0(:,2)/1e3,R0(:,3)/1e3,'LineWidth',2);
    xlabel('X(km)'); ylabel('Y(km)'); zlabel('Z(km)'); axis tight;
    hold on
    h1 = gca;
    earth_sphere(h1,'km');
    
    simula = input('Deseja simular novamente? (1 = sim) (0 = não): ');

end
    
end











