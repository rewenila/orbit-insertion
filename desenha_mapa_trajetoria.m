function desenha_mapa_trajetoria(or,traj)
% Função para desenhar um mapa em projeção, junto de uma trajetória
% fornecida
% Entradas
% or = [delta0,lamb0,h0]: origem do mapa: latitude [graus], longitude
% [graus], altitude [m]
% traj = [lat,long]: matriz onde a primeira coluna contém as latitudes e a
% segunda as longitudes de uma trajetória [graus]

ax = worldmap('World');
setm(ax,'Origin',or);
land = shaperead('landareas.shp','UseGeoCoords',true);
geoshow(land,'FaceColor',[0.15 0.8 0.15]);
cities = shaperead('worldcities','UseGeoCoords',true);
geoshow(cities,'Marker','.','Color','red');
rivers = shaperead('worldrivers','UseGeoCoords',true);
geoshow(rivers,'Color','blue');
lakes = shaperead('worldlakes','UseGeoCoords',true);
geoshow(lakes,'FaceColor','blue');
plotm(traj(1,1),traj(1,2),'*');
plotm(traj(:,1),traj(:,2));

end