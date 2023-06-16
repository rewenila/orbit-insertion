function long_c = long_ECEF2ECI(t,long,we,tg)
% Fun��o para calcular a longitude celeste a partir da longitude fixa ao
% planeta
% Entradas
% t[s]: tempo no qual se deseja saber a longitude celeste
% long [rad]: longitude relativa ao referencial fixo no planeta
% we [rad/s]: velocidade de rota��o do planeta
% tg [s]: tempo no qual o meridiano de refer�ncia tem longitude celeste
% nula
% Sa�da
% long_c [rad]: longitude celeste no tempo t

%% C�lculos
long_c = long + we*(t - tg);

end