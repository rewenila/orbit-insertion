function [D,fy,L] = aerodinamica_N_estagios(t,V,h,M,Kn,T,rho)
% Modelo de arrasto conforme a refer�ncia:
% TEWARI, A. Atmospheric and Space Flight Dynamics: Modelling and
% simulation with MATLAB and Simulink. Boston: Birkhauser, 2007.
% Exemplo 12.6
% Vale para foguete de 1, 2 ou 3 est�gios, com carga �til
% Assume-se o mesmo coeficiente de arrasto para o foguete todo e seus
% est�gios. A magnitude do arrasto � alterada em fun��o da �rea de
% refer�ncia de cada est�gio.
global ts Sr fc

%% Coeficiente de arrasto em fun��o do n�mero de Mach e de Knuden
CD = modelo_aerodinamico(V,h,M,Kn,T);
% Fator de corre��o do arrasto a partir de dados de t�nel de vento
CD = fc*CD;

%% A �rea de refer�ncia depende do est�gio atual
% N�mero de est�gios
N = length(ts);
switch N
    case 1
        S = area1estagio(t,ts,Sr);
    case 2
        S = area2estagios(t,ts,Sr);
    case 3
        S = area3estagios(t,ts,Sr);
end

%% For�as
D = 0.5*rho*V^2*S*CD;
fy = 0;
L = 0;

end 

%% �rea de refer�ncia do foguete de 1 est�gio
function S = area1estagio(t,ts,Sr)
    if t <= ts(1)
        % foguete e carga �til
        S = Sr(1);
    else
        % carga �til
        S = Sr(2);
    end
end

%% �rea de refer�ncia do foguete de 2 est�gios
function S = area2estagios(t,ts,Sr)
    if t <= ts(1)
        % todos os est�gios
        S = Sr(1);
    elseif t <= ts(2)
        % segundo est�gio e carga �til
        S = Sr(2); 
    else
        % carga �til
        S = Sr(3);
    end
end

%% �rea de refer�ncia do foguete de 3 est�gios
function S = area3estagios(t,ts,Sr)
    if t <= ts(1)
        % todos os est�gios
        S = Sr(1);
    elseif t <= ts(2)
        % segundo est�gio
        S = Sr(2); 
    elseif t <= ts(3)
        % terceiro est�gio e carga �til
        S = Sr(3); 
    else
        % carga �til
        S = Sr(4);
    end
end