function long_c = long_ECEF2ECI(t,long,we,tg)
% Função para calcular a longitude celeste a partir da longitude fixa ao
% planeta
% Entradas
% t[s]: tempo no qual se deseja saber a longitude celeste
% long [rad]: longitude relativa ao referencial fixo no planeta
% we [rad/s]: velocidade de rotação do planeta
% tg [s]: tempo no qual o meridiano de referência tem longitude celeste
% nula
% Saída
% long_c [rad]: longitude celeste no tempo t

%% Cálculos
long_c = long + we*(t - tg);

end