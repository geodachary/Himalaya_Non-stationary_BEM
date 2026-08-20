
function vl = plate_vel(LAT, LON, latp, lonp, omega) %

% function to compute local plate velocity vector given the position and the Euler pole position and angular velocity.
% INPUT:
% lat, lon - position of the point where you want velocity 
% latp, lonp, omega - Euler pole position and velocity
% latitudes and longitudes are in degrees omega is in degrees per million years
% OUTPUT:
% vl = [vn, ve, vd]' velocity in north, east, and down directions (referred to point p) in mm/yr.



Re = 6370.8e6; % Radius of Earth in millimeters
%




% convert from degrees to radians (original values are unchanged) %
deg_to_rad = pi / 180;
LAT = LAT *deg_to_rad;
LON = LON *deg_to_rad;
latp = latp*deg_to_rad;
lonp = lonp*deg_to_rad;
%
omega = omega * 1e-06 * (pi/180); 
%
% convert to Cartesian Coordinates
%

for k=1:length(LAT)
    
    lat = LAT(k);
    lon = LON(k);
    
    P = [ cos(lat)*cos(lon), cos(lat)*sin(lon), sin(lat) ]';
    EP = [cos(latp)*cos(lonp), cos(latp)*sin(lonp),sin(latp)]' * omega; %
    VC = Re * cross(EP,P); % compute the cross product: EP x P
    %
    % Rotate to local coordinate system (Cox & Hart Box 4-2);
    %
    T = zeros(3,3);
    T(1,1) = -sin(lat)*cos(lon);
    T(1,2) = -sin(lat)*sin(lon);
    T(1,3) = cos(lat);
    %
    T(2,1) = -sin(lon);

    T(2,2) = cos(lon);
    T(3,2) = 0;
    %
    T(3,1) = -cos(lat)*cos(lon); T(3,2) = -cos(lat)*sin(lon); T(3,3) = -sin(lat);
    %
    vl(:,k) = T * VC;

end
