function Gbase = Get_Gs_baselines(el,nd,xysites,rake,baselines,Vec_unit,crossind,L)


alpha = atan2(Vec_unit(:,2),Vec_unit(:,1));
        
ss = cos(rake*pi/180);
ds = sin(rake*pi/180);

xloc = [xysites';zeros(1,size(xysites,1))];


for k=1:size(el,1)
     
         
        
        
    [U, D, S] = tridisloc3d(xloc, nd', el(k,:)', [-ss(k) -ds(k) 0]', 1, .25);

   

    Veast=U(1,:)'; Vnorth=U(2,:)';
    Vbase2 = Veast(baselines(:,2)).*Vec_unit(:,1) + Vnorth(baselines(:,2)).*Vec_unit(:,2);
    Vbase1 = Veast(baselines(:,1)).*Vec_unit(:,1) + Vnorth(baselines(:,1)).*Vec_unit(:,2);
    Vbase = Vbase1 - Vbase2;

   

  
    %integrate strainrates along line for baselines that cross faults
    %(avoid velocity discontinuity)
    for j=1:size(baselines,1)
        if crossind(j)==1
            xs = xysites(baselines(j,2),1) + cos(alpha(j))*( L(j)/10/2:L(j)/10:L(j) );
            ys = xysites(baselines(j,2),2) + sin(alpha(j))*( L(j)/10/2:L(j)/10:L(j) );

          
            [U, D, S] = tridisloc3d([[xs;ys];zeros(1,length(xs))], nd', el(k,:)', [-ss(k) -ds(k) 0]', 1, .25);

            %call triangular dislocation code
            t=el(k,:);    
            [d,Strain] = CalcTriStrains_O(xs', ys', 0*ys', nd(t,1), nd(t,2), nd(t,3), 0.25, -ss(k), 0, -ds(k));    

            Exx = sum(Strain.xx);
            Exy = sum(Strain.xy);
            Eyy = sum(Strain.yy);

            Vcenters(j) = L(j)/10*(Exx*cos(alpha(j))^2 + Exy*2*cos(alpha(j))*sin(alpha(j)) + Eyy*sin(alpha(j))^2);  %equation 1.21, Segall text
 
        else

            Vcenters(j) = 0;

        end
    end

    Vbase(crossind) = Vcenters(crossind);
    Gbase(:,k)=Vbase;


    disp(['Completed ' num2str(k/size(el,1)*100) '% of baselines Greens function calculations.'] )

            
end




