function [gc,gdel] = grav_axissimetrico(r,delta)
% Entradas
% r: dist�ncia radial ao centro de massa do planeta [m]
% delta: latitude centrada no planeta [rad]
% Sa�das
% gc: componente centr�peta da gravidade [m/s^2]
% gd: componente latitudinal (norte) da gravidade [m/s^2]

%% Vari�veis
global Re mut J2 J3 J4
% Colatitude [rad]
phi = pi/2 - delta;

% Componente radial da gravidade de corpo axissim�trico
gr = -(mut)./r.^2 - (1./r).*mut*(-((J2.*Re.^2.*(-1 + 3.*cos(phi).^2))...
    ./r.^3) - (1.5.*J3.*Re.^3.*(-3.*cos(phi) + 5.*cos(phi).^3))./r.^4 - ...
    (J4.*Re.^4.*(3 - 30.*cos(phi).^2 + 35.*cos(phi).^4))./(2.*r.^5)) + ...
    (mut/r.^2).*((0.5.*J2.*Re.^2.*(-1 + 3.*cos(phi).^2))./r.^2 + (0.5...
    .*J3.*Re.^3.*(-3.*cos(phi) + 5.*cos(phi).^3))./r.^3 + (J4.*Re.^4.*...
    (3 - 30.*cos(phi).^2 + 35.*cos(phi).^4))./(8.*r.^4));

% Componente sul da gravidade de corpo axissim�trico
gphi = -(1./(r.^2)).*mut*(-((3.*J2.*Re.^2.*cos(phi).*sin(phi))./r.^2) + ...
    (0.5.*J3.*Re.^3.*(3.*sin(phi) - 15.*cos(phi).^2.*sin(phi)))./r.^3 + ...
    (J4.*Re.^4.*(60.*cos(phi).*sin(phi) - 140.*cos(phi).^3.*sin(phi)))...
    ./(8.*r.^4));

% Componentes nas dire��es centr�peta e latitudinal (norte)
gc = -gr;
gdel = -gphi;

end