function Xp = dinamica_foguete(t,X)
% Fun��o para a din�mica de transla��o de um foguete com respeito ao
% referencial PCPF
% Sistema de refer�ncia: aerodin�mico
% Sistema de coordenadas: esf�rico
% Entradas
% t: tempo [s]
% X: vetor de estado
% X(1) = V [m/s]: m�dulo do vetor velocidade relativa com respeito ao
% planeta girante
% X(2) = A [rad]: �ngulo de azimute do vetor velocidade relativa com
% respeito ao eixo z (que aponta para o norte) do sistema uen
% X(3) = phi [rad]: �ngulo de evela��o do vetor velocidade relativa com
% respeito ao horizonte local (plano yz do referencial uen)
% X(4) = r [m]: dist�ncia radial at� o centro do planeta
% X(5) = delta [rad]: latitude com respeito ao plano equatorial do planeta
% X96) = lon (rad): longitude planet�ria
% Sa�da
% Xp: derivada do vetor de estado X

%% Entrada de constantes por vari�veis globais
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

%% Fun��o para c�lculo da massa e tra��o em fun��o do tempo
% Depende da taxa de queima de propelente, tempo de queima de propelente,
% impulso espec�fico e separa��o de est�gios
[ft,m] = propulsao_N_estagios(t);

%% Os �ngulos de apontamento da tubeira s�o nulos
epsl = 0; mu = 0;

%% Fun��o para c�lculo do modelo atmosf�rico
h = r - Re; % altitude
[T,~,rho,~,M,~,~,Kn,~,~] = atm_padrao(h,V,lc,dT);

%% Fun��o para c�lculo do modelo aerodin�mico
% Depende da altitude e velocidade
[D,fy,L] = aerodinamica_N_estagios(t,V,h,M,Kn,T,rho);

%% C�lculo da gravidade
% Fun��o para c�lculo do modelo gravitacional
[gc,gd] = grav_axissimetrico(r,delta);

%% Equa��es de cinem�tica de transla��o
rp = V*sin(phi);
deltap = (V/r)*cos(phi)*cos(A);
lonp = (V*cos(phi)*sin(A))/(r*cos(delta));

%% Equa��es de din�mica de transla��o
Vp = (1/m)*(ft*cos(epsl)*cos(mu) - D - m*gc*sin(phi) + m*gd*cos(phi)*cos(A)...
    - m*we^2*r*cos(delta)*(cos(phi)*cos(A)*sin(delta) - sin(phi)*cos(delta)));

Ap = (1/(m*V*(cos(phi))))*(m*(V^2/r)*cos(phi)^2*sin(A)*tan(delta)...
    + ft*sin(mu) + fy - m*gd*sin(A) + m*we^2*r*sin(A)*sin(delta)*cos(delta)...
    - 2*m*we*V*(sin(phi)*cos(A)*cos(delta) - cos(phi)*sin(delta)));

phip = (1/(m*V))*(m*(V^2/r)*cos(phi) + ft*sin(epsl)*cos(mu) + L - m*gc*cos(phi)...
    - m*gd*sin(phi)*cos(A) + m*we^2*r*cos(delta)*sin(sin(phi)*cos(A)*sin(delta)...
    + cos(phi)*cos(delta)) + 2*m*we*V*sin(A)*cos(delta));

%% Satura��o da altitude
if h < 0 % altitude negativa n�o � permitida
    % mant�m as derivadas nulas
    rp = 0; deltap = 0; lonp = 0; Vp = 0; Ap = 0; phip = 0;
end

%% Modela o trilho de lan�amento
H = h - h0; % altura
if ((H <= l_trilho) && (t <= 10)) % verifica se a altura � menor que l_trilho nos primeiros segundos da simula��o
    Ap = 0; phip = 0; % anula as derivadas dos �ngulos de orienta��o da velocidade
end

%% Derivadas do vetor de estado
Xp = [Vp;Ap;phip;rp;deltap;lonp];

end