function [D,fy,L] = aerodinamica_N_estagios(t,V,h,M,Kn,T,rho)
% Modelo de arrasto conforme a referência:
% TEWARI, A. Atmospheric and Space Flight Dynamics: Modelling and
% simulation with MATLAB and Simulink. Boston: Birkhauser, 2007.
% Exemplo 12.6
% Vale para foguete de 1, 2 ou 3 estágios, com carga útil
% Assume-se o mesmo coeficiente de arrasto para o foguete todo e seus
% estágios. A magnitude do arrasto é alterada em função da área de
% referência de cada estágio.
global ts Sr fc

%% Coeficiente de arrasto em função do número de Mach e de Knuden
CD = modelo_aerodinamico(V,h,M,Kn,T);
% Fator de correção do arrasto a partir de dados de túnel de vento
CD = fc*CD;

%% A área de referência depende do estágio atual
% Número de estágios
N = length(ts);
switch N
    case 1
        S = area1estagio(t,ts,Sr);
    case 2
        S = area2estagios(t,ts,Sr);
    case 3
        S = area3estagios(t,ts,Sr);
end

%% Forças
D = 0.5*rho*V^2*S*CD;
fy = 0;
L = 0;

end 

%% Área de referência do foguete de 1 estágio
function S = area1estagio(t,ts,Sr)
    if t <= ts(1)
        % foguete e carga útil
        S = Sr(1);
    else
        % carga útil
        S = Sr(2);
    end
end

%% Área de referência do foguete de 2 estágios
function S = area2estagios(t,ts,Sr)
    if t <= ts(1)
        % todos os estágios
        S = Sr(1);
    elseif t <= ts(2)
        % segundo estágio e carga útil
        S = Sr(2); 
    else
        % carga útil
        S = Sr(3);
    end
end

%% Área de referência do foguete de 3 estágios
function S = area3estagios(t,ts,Sr)
    if t <= ts(1)
        % todos os estágios
        S = Sr(1);
    elseif t <= ts(2)
        % segundo estágio
        S = Sr(2); 
    elseif t <= ts(3)
        % terceiro estágio e carga útil
        S = Sr(3); 
    else
        % carga útil
        S = Sr(4);
    end
end