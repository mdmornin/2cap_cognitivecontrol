% Offset X-Coordinate Based on Kcoords
shank1 = kcoords == 1;
shank2 = kcoords == 2;




xcoords(shank1) = xcoords(shank1) + 0;
xcoords(shank2) = xcoords(shank2) + 1000;


