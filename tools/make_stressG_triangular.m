function stressG = make_stressG_triangular(tri,nd,rake)



patch_stuff = make_triangular_patch_stuff(tri,nd);


centroids = patch_stuff.centroids_faces;
NormalFaces = patch_stuff.normal_faces;
StrikeVecFaces = patch_stuff.strikevec_faces;


%create rake vector as a rotation of strike vector about normal
   
R = zeros(3,3);

for k=1:size(NormalFaces,1)
    n = NormalFaces(k,:);  %negative so normal points up
    a = rake(k)*pi/180;  
    R(1,:) = [cos(a)+n(1)^2*(1-cos(a))          n(1)*n(2)*(1-cos(a))-n(3)*sin(a)  n(1)*n(3)*(1-cos(a))+n(2)*sin(a)];
    R(2,:) = [n(1)*n(2)*(1-cos(a))+n(3)*sin(a)  cos(a)+n(2)^2*(1-cos(a))          n(2)*n(3)*(1-cos(a))-n(1)*sin(a)];
    R(3,:) = [n(1)*n(3)*(1-cos(a))-n(2)*sin(a)  n(2)*n(3)*(1-cos(a))+n(1)*sin(a)  cos(a)+n(3)^2*(1-cos(a))];

    rakevec(k,:) = (R*StrikeVecFaces(k,:)')';

end


    
    
    

stressG = zeros(size(tri,1),size(tri,1));


ss = cos(rake*pi/180);
ds = sin(rake*pi/180);

mu = 3*10^10;
nd = 1000*nd;
centroids = centroids*1000;


for j=1:size(tri,1)


    temp1{1}=nd;
    temp2{1}=tri(j,:);    



    %[U, D, S] = tridisloc3d(centroids', temp1{1}', temp2{1}', [ss(j) -ds(j) 0]', mu, .25);  %negative indicates reverse slip

    [d,Strain] = CalcTriStrains_O(centroids(:,1), centroids(:,2), centroids(:,3), temp1{1}(temp2{1},1), temp1{1}(temp2{1},2), temp1{1}(temp2{1},3), 0.25, ss(j), 0, -ds(j));    
    S = StrainToStress(Strain, mu, mu);
     

   
    %calculate tractions
    T1=S.xx'.*NormalFaces(:,1)+S.xy'.*NormalFaces(:,2)+S.xz'.*NormalFaces(:,3);
    T2=S.xy'.*NormalFaces(:,1)+S.yy'.*NormalFaces(:,2)+S.yz'.*NormalFaces(:,3);
    T3=S.xz'.*NormalFaces(:,1)+S.yz'.*NormalFaces(:,2)+S.zz'.*NormalFaces(:,3);
   
    %T1=S(1,:)'.*NormalFaces(:,1)+S(2,:)'.*NormalFaces(:,2)+S(3,:)'.*NormalFaces(:,3);
    %T2=S(2,:)'.*NormalFaces(:,1)+S(4,:)'.*NormalFaces(:,2)+S(5,:)'.*NormalFaces(:,3);
    %T3=S(3,:)'.*NormalFaces(:,1)+S(5,:)'.*NormalFaces(:,2)+S(6,:)'.*NormalFaces(:,3);
   
    
    %component of traction in rake direction
    GTv = T1.*rakevec(:,1) + T2.*rakevec(:,2) + T3.*rakevec(:,3);
     
  
    stressG(:,j) = GTv;   %negative sign so that stress drops with slip
   


end %j



