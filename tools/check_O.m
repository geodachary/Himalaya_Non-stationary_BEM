function [sx,sy,sz] = check_O(sx, sy, sz, X, Y, Z,small)
%This function is to avoid the legs of the triangular dislocation.
%input in global coordinate system
%U: displacement 3xN matrix in global coordinate system

%Transform from global coordinate system to local coordinate system 
sy = -sy;
sz = -sz;
Y = -Y;
Z = -Z;

%Check if any point is above free surface 
%Check if any point is above free surface 
clear index
index = Z < 0;
if sum(index) > 0
    disp(' input of Z must be all negative')
    disp(' There is a pause in check_O');
    pause
end

clear index
index = sz < 0;
if sum(index) >0
    disp(' input of sz must be all negative')
    disp(' There is a pause in check_O');
    pause
end

%Check if the point is on any of leg plane.
%leg 1
plane1 =(Y(1)-Y(2))*sx+(-X(1)+X(2))*sy + (Y(1)-Y(2))*X(1)+(-X(1)+X(2))*Y(1);
index1 = plane1 == 0;
if sum(index1)> 0
    disp('hit leg')
sx(index1)= sx(index1) + small;
sy(index1)= sy(index1) + small;
end
%leg 2
plane2 =(Y(2)-Y(3))*sx+(-X(2)+X(3))*sy + (Y(2)-Y(3))*X(2)+(-X(2)+X(3))*Y(2);
index2 = plane2 == 0;
if sum(index2)> 0
sx(index2)= sx(index2) + small;
sy(index2)= sy(index2) + small;
end
%leg 3
plane3 =(Y(3)-Y(1))*sx+(-X(3)+X(1))*sy + (Y(3)-Y(1))*X(3)+(-X(3)+X(1))*Y(3);
index3 = plane3 == 0;
if sum(index3)> 0
sx(index3)= sx(index3) + small;
sy(index3)= sy(index3) + small;
end
%Coordinate transformation
  sy =-sy;
  sz =-sz;
  

