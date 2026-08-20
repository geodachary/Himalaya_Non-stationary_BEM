function [sx, sy, sz, index_sp] = CalcTriDisps_Check(sx, sy, sz, X, Y, Z, pr, ss, ts, ds, small);
%This function is to make sure that the observed positions are not within
%the range of any singularity of small so that the strain and stress can be
%accurate enough.

%Author: Wen-Jeng Huang
%date: August 25th, 2007
%latest revision: August 26th, 2007

[sx,sy,sz] = check_O(sx, sy, sz, X, Y, Z, 2*small);
%input in global coordinate system
%U: displacement 3xN matrix in global coordinate system


%transform from global coordinate system to local coordinate system 
sy = -sy;
sz = -sz;
Y = -Y;
Z = -Z;
%change slip sense to fit Okada's convetion
ts = -ts;
ds = -ds;

%Avoid some special singularities If there is one, shit 0.01 in x direction
[u,index_sp] = CalcTriDisps_Ch(sx, sy, sz, X, Y, Z, pr, ss, ts, ds);
     if sum(index_sp)>0
       sx(logical(index_sp))= sx(logical(index_sp))+ 0.01;
       %disp('Shifting occurs')
     end


%Check if displacements in any position are real numbers.
[u] = CalcTriDisps_Ch(sx, sy, sz, X, Y, Z, pr, ss, ts, ds);

U =[u.x,u.y,u.z]';
clear u index
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

if sum(index)>0 % if sum(index)>0

    index = logical(index);
    if sum(index) > 0
        %disp('Shifting occurs')
         %disp('nan or complex number occurs')
    v1 =(Y(1)-Y(2))*(Z(3)-Z(2))-(Y(3)-Y(2))*(Z(1)-Z(2));
    v2 =(Z(1)-Z(2))*(X(3)-X(2))-(Z(3)-Z(2))*(X(1)-X(2));
    v3 =(X(1)-X(2))*(Y(3)-Y(2))-(X(3)-X(2))*(Y(1)-Y(2));
    vector =[v1 v2 v3];
    vecLeg=[X(1) Y(1) Z(1)];
    vec = 1.1*small*cross(vector,vecLeg)/norm(cross(vector,vecLeg));
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
 
end
 
 sy = -sy;
 sz = -sz;