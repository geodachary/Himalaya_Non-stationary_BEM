function [d,strain] = CalcTriStrains_O(x1, x2, x3, x, y, z, pr, ss, ts, ds)
%This function is to calulate strain
%It need to use with several additional functions in order to make sure
%   that the observed position outside of the range of any singularity of small(0.0001 unit).
%input in global coordinate system
%output in global coordinate sytem as well
%d: Nx9 matrix, where N is the amount of calculated points.
%strain: strain.xx, 1*N matrix;
%        strain.xy  1*N matrix;
%        strain.xz  1*N matrix;
%        strain.yy  1*N matrix;
%        strain.yz  1*N matrix;
%        strain.zz  1*N matrix;

%Author: Wen-Jeng Huang
%date: August 25th, 2007
%latest revison: August 30th, 2007

    %Not need for now
    %v1 =(y(1)-y(2))*(z(3)-z(2))-(y(3)-y(2))*(z(1)-z(2));
    %v2 =(z(1)-z(2))*(x(3)-x(2))-(z(3)-z(2))*(x(1)-x(2));
    %v3 =(x(1)-x(2))*(y(3)-y(2))-(x(3)-x(2))*(y(1)-y(2));
    %vector =[v1 v2 v3];
    %area = 1/2*norm(vector);
    %small_normalize = 0.0001*area;
 
    %small: half of length of a cubic which is used to caculate strain 
 small = 0.0001;   
%check if the point is the centroid of the triangle.
%if it is, move 1.01 small for calulating strain.
centroidx =1/3*sum(x);
centroidy =1/3*sum(y);
centroidz =1/3*sum(z);

index1 = abs(x1 - centroidx) <0.00001;
index2 = abs(x2 - centroidy) <0.00001;
index3 = abs(x3 - centroidz) <0.00001;

index =index1 + index2 + index3;

index4= index == 3;
if sum(index4) ==0
    
elseif sum(index4) == 1
    v1 =(y(1)-y(2))*(z(3)-z(2))-(y(3)-y(2))*(z(1)-z(2));
    v2 =(z(1)-z(2))*(x(3)-x(2))-(z(3)-z(2))*(x(1)-x(2));
    v3 =(x(1)-x(2))*(y(3)-y(2))-(x(3)-x(2))*(y(1)-y(2));
    vector =[v1 v2 v3];
    vector= 1.01*small*(vector/norm(vector));
    x1(index1)= centroidx +vector(1);
    x2(index2)= centroidy +vector(2);
    x3(index3)= centroidz +vector(3);
    
    %Avoid to be above free surface
        if vector(3)<0
    x1(index1)= centroidx +vector(1);
    x2(index2)= centroidy +vector(2);
    x3(index3)= centroidz +vector(3);
        else
    x1(index1)= centroidx -vector(1);
    x2(index2)= centroidy -vector(2);
    x3(index3)= centroidz -vector(3);
        end    
else
    disp('There is more than one point on the centroid')
    disp('There is a pause in CalcTriStrains(line 67)')
    pause
end


% Check if any observed points are within the range of any singularity.
%If there are, shift the points some small distance.
[x1, x2, x3,index_sp] = CalcTriDisps_Check(x1, x2, x3, x, y, z, pr, ss, ts, ds,small);

%Check if all point is below depth of -small; Make sure all calculated
%points are below free surface
index = x3 >=-small;
if sum(index)>0
    x3(index)= x3(index)-small;
end
   xcubic =zeros(6,size(x1,1));
   ycubic =zeros(6,size(x1,1));
   zcubic =zeros(6,size(x1,1));

for counter = 1 :6 %A4
      
    switch counter  % A5
        case 1
           x1n = x1  + small;
           x2n = x2;
           x3n = x3 ; 
       case 2
           x1n = x1;
           x2n = x2 + small ;
           x3n = x3 ;
       case 3
           x1n = x1; 
           x2n = x2;
           x3n = x3 +small;
       case 4
           x1n = x1 -small;
           x2n = x2 ;
           x3n = x3;
       case 5
           x1n = x1;
           x2n = x2 - small;
           x3n = x3;
       case 6 
           x1n = x1;
           x2n = x2;
           x3n = x3 - small;
   end          %A5
  
% Save time if the fault is pure dipslip or pure strike slip  
   if   ss == 0 && ts==0 && ds == 0
       u=[zeros(size(x1n,1),1);zeros(size(x1n,1),1);zeros(size(x1n,1),1)];
   else     
    [u] = CalcTriDisps_O(x1n, x2n, x3n, x, y, z, pr, ss, ts, ds);
   end
     
  xcubic(counter,:)= u(1,:);
  ycubic(counter,:)= u(2,:);
  zcubic(counter,:)= u(3,:);
  
end
  %
  Dxu=zeros(1,size(xcubic,2));
  Dxv=zeros(1,size(xcubic,2));
  Dxw=zeros(1,size(xcubic,2));
  Dyu=zeros(1,size(xcubic,2));
  Dyv=zeros(1,size(xcubic,2));
  Dyw=zeros(1,size(xcubic,2));
  Dzu=zeros(1,size(xcubic,2));
  Dzv=zeros(1,size(xcubic,2));
  Dzw=zeros(1,size(xcubic,2));
  
  Dxu =(xcubic(1,:)-xcubic(4,:))/(2*small);
  Dxv =(ycubic(1,:)-ycubic(4,:))/(2*small);
  Dxw =(zcubic(1,:)-zcubic(4,:))/(2*small);
  
  Dyu =(xcubic(2,:)-xcubic(5,:))/(2*small);
  Dyv =(ycubic(2,:)-ycubic(5,:))/(2*small);
  Dyw =(zcubic(2,:)-zcubic(5,:))/(2*small);
  
  Dzu =(xcubic(3,:)-xcubic(6,:))/(2*small);
  Dzv =(ycubic(3,:)-ycubic(6,:))/(2*small);
  Dzw =(zcubic(3,:)-zcubic(6,:))/(2*small);
  
  %displacement gradients
  d=[Dxu;Dyu;Dzu;Dxv;Dyv;Dzv;Dxw;Dyw;Dzw]';
  
  strain.xx = Dxu;
  strain.xy = 1/2*(Dyu+Dxv);
  strain.xz = 1/2*(Dzu+Dxw);
  strain.yy = Dyv;
  strain.yz = 1/2*(Dyw+Dzv);
  strain.zz = Dzw;  

 