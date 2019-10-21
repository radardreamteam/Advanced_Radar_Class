function ecef = lla2ecef(lla)
% Usage: ecef = lla2ecef(lla)
% Converts lla coodinate matrix in to earth centered earth fixed. The units of
% ecef are the same as the units of altitude provided in lla.
% lla = [latitude (degrees) longitude (degrees) altitude (arb)]
%
    [x,y,z]=geodetic2ecef(lla(:,1)*pi/180,lla(:,2)*pi/180,lla(:,3),referenceEllipsoid('wgs84'));
    ecef=[x y z];
end