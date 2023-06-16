function Xp = dinamica_foguete(t,X)
% Função para a dinâmica de translação de um foguete com respeito ao
% referencial PCPF
% Sistema de referência: aerodinâmico
% Sistema de coordenadas: esférico
% Entradas
% t: tempo [s]
% X: vetor de estado
% X(1) = V [m/s]: módulo do vetor velocidade relativa com respeito ao
% planeta girante
% X(2) = A [rad]: ângulo de azimute do vetor velocidade relativa com
% respeito ao eixo z (que aponta para o norte) do sistema uen
% X(3) = phi [rad]: ângulo de evelação do vetor velocidade relativa com
% respeito ao horizonte local (plano yz do referencial uen)
% X(4) = r [m]: distância radial até o centro do planeta
% X(5) = delta [rad]: latitude com respeito ao plano equatorial do planeta
% X96) = lon (rad): longitude planetária
% Saída
% Xp: derivada do vetor de estado X

%% Entrada de constantes por variáveis globais
global we Re lc dT h0 l_trilho

%% Vetor de estado
V = X(1);
A = X(2);
phi = X(3);
r = X(4);
delta = X(5);
%lon = X(6);
if V < 0
    V = 0; % evita velocidade negativa
end

%% Função para cálculo da massa e tração em função do tempo
% Depende da taxa de queima de propelente, tempo de queima de propelente,
% impulso específico e separação de estágios
[ft,m] = propulsao_N_estagios(t);

%% Os ângulos de apontamento da tubeira são nulos
epsl = 0; mu = 0;

%% Função para cálculo do modelo atmosférico
h = r - Re; % altitude
[T,~,rho,~,M,~,~,Kn,~,~] = atm_padrao(h,V,lc,dT);

%% Função para cálculo do modelo aerodinâmico
% Depende da altitude e velocidade
[D,fy,L] = aerodinamica_N_estagios(t,V,h,M,Kn,T,rho);

%% Cálculo da gravidade
% Função para cálculo do modelo gravitacional
[gc,gd] = grav_axissimetrico(r,delta);

%% Equações de cinemática de translação
rp = V*sin(phi);
deltap = (V/r)*cos(phi)*cos(A);
lonp = (V*cos(phi)*sin(A))/(r*cos(delta));

%% Equações de dinâmica de translação
Vp = (1/m)*(ft*cos(epsl)*cos(mu) - D - m*gc*sin(phi) + m*gd*cos(phi)*cos(A)...
    - m*we^2*r*cos(delta)*(cos(phi)*cos(A)*sin(delta) - sin(phi)*cos(delta)));

Ap = (1/(m*V*(cos(phi))))*(m*(V^2/r)*cos(phi)^2*sin(A)*tan(delta)...
    + ft*sin(mu) + fy - m*gd*sin(A) + m*we^2*r*sin(A)*sin(delta)*cos(delta)...
    - 2*m*we*V*(sin(phi)*cos(A)*cos(delta) - cos(phi)*sin(delta)));

phip = (1/(m*V))*(m*(V^2/r)*cos(phi) + ft*sin(epsl)*cos(mu) + L - m*gc*cos(phi)...
    - m*gd*sin(phi)*cos(A) + m*we^2*r*cos(delta)*sin(sin(phi)*cos(A)*sin(delta)...
    + cos(phi)*cos(delta)) + 2*m*we*V*sin(A)*cos(delta));

%% Saturação da altitude
if h < 0 % altitude negativa não é permitida
    % mantém as derivadas nulas
    rp = 0; deltap = 0; lonp = 0; Vp = 0; Ap = 0; phip = 0;
end

%% Modela o trilho de lançamento
H = h - h0; % altura
if ((H <= l_trilho) && (t <= 10)) % verifica se a altura é menor que l_trilho nos primeiros segundos da simulação
    Ap = 0; phip = 0; % anula as derivadas dos ângulos de orientação da velocidade
end

%% Derivadas do vetor de estado
Xp = [Vp;Ap;phip;rp;deltap;lonp];

end