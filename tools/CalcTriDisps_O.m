function [U] = CalcTriDisps_O(sx, sy, sz, X, Y, Z, pr, ss, ts, ds)
%This function is to calculate displacements
%It need to be used with another function, CalctriDisps_Ch

%Author: Wen-Jeng Huang
%date:August 25th, 2007
%latest revision: Auguest 30th, 2007

%sx, sy, sz are observed positions in a global coordinate system. Poitive direction of axese are
%east(sx), north(sy) and up(sz).
%X, Y, Z are vertices of triangular dislocation in global coordinate.
%U is a matrix of displacement in global coordinate system.

 small = 0.0001;%For shifting the position in order to avoid singularity
% Check if any position on the legs
% If it is, shift 0.0002
[sx,sy,sz] = check_O(sx, sy, sz, X, Y, Z,2*small);

%Input in global coordinate system
%U: displacement 3xN matrix in global coordinate system
%Transform from global coordinate system to local coordinate system 
sy = -sy;
sz = -sz;
Y = -Y;
Z = -Z;
%Change slip sense to fit Okada's convetion
ts = -ts;
ds = -ds;

%To calculate displacement
%All input and output are in local coordinate system;
%sx(positive) in east, sy(positve)in south, sz(positive) in north
[u,index] = CalcTriDisps_Ch(sx, sy, sz, X, Y, Z, pr, ss, ts, ds);
%index is a indicator for finding singularity      
           indextemp = index ==0;
               if sum(index)>0
       sx(logical(index))= sx(logical(index)) + 2*small;
       Newsx =sx; Newsy =sy; Newsz = sz;
       Newsx(indextemp)=[];Newsy(indextemp)=[];Newsz(indextemp)=[];
       [uu] = CalcTriDisps_Ch(Newsx, Newsy, Newsz, X, Y, Z, pr, ss, ts, ds);
       u.x(logical(index))=uu.x;
       u.y(logical(index))=uu.y;
       u.z(logical(index))=uu.z;
               end    
U =[u.x,u.y,u.z]';

%check if there is nan or complex number.
%If there is, shift the point small in the directon of the normal vector of
%the dislocation
clear u index indextemp
for loop = 1: size(U,2)
    index1(loop) = isnan(U(1,loop));
    index2(loop) = isnan(U(2,loop));
    index3(loop) = isnan(U(3,loop));

    if isreal(U(1,loop)) == 0
    index1(loop) = 1;
    end
    if isreal(U(2,loop)) == 0
    index2(loop) = 1;
    end
    if isreal(U(2,loop)) == 0
    index3(loop) = 1;
    end
end

index = index1 + index2 + index3;
if sum(index)>0
    index = logical(index);

    if sum(index) > 0
    v1 =(Y(1)-Y(2))*(Z(3)-Z(2))-(Y(3)-Y(2))*(Z(1)-Z(2));
    v2 =(Z(1)-Z(2))*(X(3)-X(2))-(Z(3)-Z(2))*(X(1)-X(2));
    v3 =(X(1)-X(2))*(Y(3)-Y(2))-(X(3)-X(2))*(Y(1)-Y(2));
    vector =[v1 v2 v3];
    vecLeg=[X(1) Y(1) Z(1)];
    vec = small*cross(vector,vecLeg)/norm(cross(vector,vecLeg));
    
    %To avoid the sz might be above free surface
     if vec(3)> 0
    sx(index)= sx(index) + vec(1);
    sy(index)= sy(index) + vec(2);
    sz(index)= sz(index) + vec(3);
     else
    sx(index)= sx(index) - vec(1);
    sy(index)= sy(index) - vec(2);
    sz(index)= sz(index) - vec(3);
     end 
    end
        
    [u] = CalcTriDisps_Ch(sx, sy, sz, X, Y, Z, pr, ss, ts, ds);
    Ut =[u.x,u.y,u.z]';
    index = logical([index;index;index]);
    indextemp = index ==0;
    Ut(indextemp)=[];
    U(index)= Ut;   
end
    
%check if 
clear u index index1 index2 index3
for loop = 1: size(U,2)
    index1(loop) = isnan(U(1,loop));
    index2(loop) = isnan(U(2,loop));
    index3(loop) = isnan(U(3,loop));

    if isreal(U(1,loop)) == 0
    index1(loop) = 1;
    end
    if isreal(U(2,loop)) == 0
    index2(loop) = 1;
    end
    if isreal(U(2,loop)) == 0
    index3(loop) = 1;
    end
end

index = index1 + index2 + index3;
 if sum(index) >0
 disp('There is still something wrong')
 disp('There is a pause in CalcTriDisps_O (line 108).')
 pause  
 end
 
 %Return to Global coordinate system
 U =U'*[1,0,0;0,-1,0;0,0,-1];
 U =U';
