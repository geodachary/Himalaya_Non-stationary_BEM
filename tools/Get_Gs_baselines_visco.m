function Gbase = Get_Gs_baselines_visco(Ge, Gn, xysites, baselines, Vec_unit, crossind, L, el, nd, rake)

    alpha = atan2(Vec_unit(:,2), Vec_unit(:,1));
    ss = cos(rake * pi / 180);
    ds = sin(rake * pi / 180);
    
    Gbase = zeros(size(baselines,1), size(el,1));
    
    for k = 1:size(el,1)
        Veast = Ge(:,k);
        Vnorth = Gn(:,k);
        
        Vbase1 = Veast(baselines(:,1)) .* Vec_unit(:,1) + Vnorth(baselines(:,1)) .* Vec_unit(:,2);
        Vbase2 = Veast(baselines(:,2)) .* Vec_unit(:,1) + Vnorth(baselines(:,2)) .* Vec_unit(:,2);
        Vbase = Vbase1 - Vbase2;
        
        % Strain rate integration for fault-crossing baselines
        Vcenters = zeros(size(baselines,1), 1);
        for j = 1:size(baselines,1)
            if crossind(j) == 1
                xs = xysites(baselines(j,2),1) + cos(alpha(j)) * (L(j)/10/2 : L(j)/10 : L(j));
                ys = xysites(baselines(j,2),2) + sin(alpha(j)) * (L(j)/10/2 : L(j)/10 : L(j));
                
                t = el(k,:);
                [d, Strain] = CalcTriStrains_O(xs', ys', 0*ys', nd(t,1), nd(t,2), nd(t,3), 0.25, -ss(k), 0, -ds(k));
                
                Exx = sum(Strain.xx);
                Exy = sum(Strain.xy);
                Eyy = sum(Strain.yy);
                
                Vcenters(j) = L(j)/10 * (Exx * cos(alpha(j))^2 + Exy * 2 * cos(alpha(j)) * sin(alpha(j)) + Eyy * sin(alpha(j))^2);
            end
        end
        
        Vbase(crossind) = Vcenters(crossind);
        Gbase(:,k) = Vbase;
        
        disp(['Completed ' num2str(k/size(el,1)*100) '% of baselines Greens function calculations.'] )
    end
end



