function par_orb = det_orb(t0,rc0,vc0,mu)
% Função para determinar parâmetros orbitais a partir de uma observação de
% posição e outra de velocidade, sendo as mesmas tomadas em relação ao
% primário em um problema de dois corpos e escritas no referencial
% celestial
% Entradas
% t0: tempo em que a observação foi feita (s)
% rc0: vetor posição relativa de m2 com respeito a m1 escrito no
% referencial celeste (m ou km)
% vc0: vetor velocidade relativa escrito no referencial celeste (m/s ou
% km/s, unidades coerentes com o vetor posição); deve ser tomado no mesmo
% instante da medida de posição
% mu: parâmetro gravitacional padrão do corpo m1 (m^3/s^2 ou km^3/s^2,
% unidades coerentes com a posição e velocidade)
% Saídas
% par_orb: vetor de parâmetros orbitais
% a = par_orb(1): semi eixo maior da órbita (m ou km, depende das unidades
% de entrada
% e = par_orb(2): excentricidade da órbita (adimensional)
% tau = par_orb(3): tempo de periastro (s)
% OMEGA = par_orb(4): longitude celeste do nodo ascendente (rad)
% i = par_orb(5): inclinação (rad)
% omega = par_orb(6): argumento de periastro (rad)

%% Cálculos
% Distância radial no instante observado
r0 = norm(rc0);
% Vetor quantidade de movimento angular específica no referencial celeste
hc = cross(rc0,vc0);
% Vetor escentricidade no referencial celeste
ec = cross(vc0,hc)/mu - rc0/r0;
% Excentricidade da órbita
e = norm(ec);
% Módulo do vetor hc
h = norm(hc);
% Parâmetro da órbita
p = h^2/mu;
% Semi eixo maior
a = p/(1 - e^2);
% Vetor parâmetro no referencial celeste
pc = p*cross(hc,ec)/(h*e);
% Anomalia verdadeira no instante da observação
costheta = (p - r0)/(e*r0);
sintheta = dot(rc0,pc)/(r0*p);
theta0 = atan2(sintheta,costheta);

% O tempo de periastro depende do tipo de órbita
if (0 <= e) && (e < 1)
    tipo = 'e'; % órbita elíptica
elseif e == 1
    tipo = 'p'; % órbita parabólica
else
    tipo = 'h'; % órbita hiperbólica
end

% Tempo de periastro
if tipo == 'e'
    % Movimento médio
    n = sqrt(mu/a^3);
    % Anomalia excêntrica
    E0 = 2*atan(sqrt((1 - e)/(1 + e))*tan(theta0/2));
    % Tempo de periastro
    tau = t0 - (E0 - e*sin(E0))/n;
elseif tipo == 'p'
    tau = t0 - ((tan(theta0/2))^3 + 3*tan(theta0/2))/(mu/p^3)^(1/6);
else
    % Anomalia hiperbólica
    H0 = 2*atanh(sqrt((e - 1)/(1 + e))*tan(theta0/2));
    % Tempo de periastro
    tau = t0 - (e*sinh(H0) - H0)/((e^2 - 1)^(3/2)*sqrt(mu/p^3));
end

% Linha dos nodos
% Vetor unitário ao longo d alinha dos nodos (no sistema celeste)
Kc = [0;0;1];
nc = cross(Kc,hc)/norm(cross(Kc,hc));
% Longitude celeste do nodo ascendente
OMEGA = atan2(nc(2),nc(1));
% Inclinação
if abs(OMEGA) > 1e-3 % evita divisão por zero ou por um número perto de zero
    i = atan2(hc(1)/sin(OMEGA),hc(3));
else
    i = atan2(-hc(2)/cos(OMEGA),hc(3));
end    
% Vetor unitário ao longo do vetor excentricidade (no referencial celeste)
ie = ec/e;
% Vetor unitário ao longo do vetor h (no sistema celeste)
ih = hc/h;
% Argumento de periastro
cosomega = dot(ie,nc);
sinomega = dot(ih,cross(nc,ie));
omega = atan2(sinomega,cosomega);

%% Vetor parâmetros de saída
par_orb(1) = a;
par_orb(2) = e;
par_orb(3) = tau;
par_orb(4) = OMEGA;
par_orb(5) = i;
par_orb(6) = omega;
    
end 