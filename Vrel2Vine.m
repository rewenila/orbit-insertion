function [vi,phii,Ai] = Vrel2Vine(vr,phir,Ar,we,r,dt)
% Fun��o que converte a velocidade, eleva��o e azimute da velocidade
% relativa para a velocidade inercial
% Entradas
% vr [m/s]: velocidade relativa
% phir [rad]: inclina��o da velocidade relativa
% Ar [rad]: azimute da velocidade relativa
% we [rad/s]: velocidade de rota��o do referencial girante
% r[m]: dist�ncia radial at� a origem do referencial inercial
% dt [rad]: latitude
% Sa�das
% v [m/s]: magnitude da velocidade com respeito ao referencial inercial
% phi [rad]: �ngulo de eleva��o da velocidade inercial (�ngulo de
% trajet�ria)
% A [rad]: �ngulo de azimute da velocidade inercial

%% C�lculos
Ai = atan2(vr*cos(phir)*sin(Ar) + we*r*cos(dt),vr*cos(phir)*cos(Ar));
if Ai < 0
    Ai = Ai + 2*pi;
end
vi = sqrt(vr^2 + 2*vr*cos(phir)*sin(Ar)*r*we*cos(dt) + r^2*we^2*cos(dt)^2);
phii = atan2(sin(phir)*cos(Ai),cos(phir)*cos(Ar));

end