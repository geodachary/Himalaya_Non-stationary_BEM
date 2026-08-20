
function [Vbase,Sigbase,BaseEnds,L,baselines,Vec_unit,crossind] = make_baseline_rate_changes(xy,Veast,Vnorth,Sigeast,Signorth,PatchEnds)



if nargin<6
    PatchEnds = [];
end


%plot data and triangles
 tri = delaunay(xy(:,1),xy(:,2));
 %remove triangles with small angles

for k=1:size(tri,1)
    v1 = [xy(tri(k,3),1) xy(tri(k,3),2)]-[xy(tri(k,1),1) xy(tri(k,1),2)]; 
    v2 = [xy(tri(k,2),1) xy(tri(k,2),2)]-[xy(tri(k,1),1) xy(tri(k,1),2)];
    a1 = 180/pi*acos( dot(v1/norm(v1),v2/norm(v2)) );
    v1 = [xy(tri(k,1),1) xy(tri(k,1),2)]-[xy(tri(k,2),1) xy(tri(k,2),2)]; 
    v2 = [xy(tri(k,3),1) xy(tri(k,3),2)]-[xy(tri(k,2),1) xy(tri(k,2),2)];
    a2 = 180/pi*acos( dot(v1/norm(v1), v2/norm(v2)) );
    v1 = [xy(tri(k,1),1) xy(tri(k,1),2)]-[xy(tri(k,3),1) xy(tri(k,3),2)]; 
    v2 = [xy(tri(k,2),1) xy(tri(k,2),2)]-[xy(tri(k,3),1) xy(tri(k,3),2)];
    a3 = 180/pi*acos( dot(v1/norm(v1),v2/norm(v2)) );
    
    if abs(a1)<10 | abs(a2)<10 | abs(a3)<10
        index(k)=1;
    else
        index(k)=0;
    end

end
index=logical(index);
tri(index,:)=[];


%form baselines
tri_num = (1:size(tri,1))';
tri_num = repmat(tri_num,3,1);
baselines = [ tri(:,1:2) ; tri(:,2:3) ; [tri(:,1) tri(:,3) ]];  %keep track of triangle number
baselines = sort(baselines,2);

%baseline endpoints
BaseEnds = [xy(baselines(:,1),:) xy(baselines(:,2),:)];

%identify baslines that cross patches (will use model strainrate value at center of baseline)



if ~isempty(PatchEnds)
    for k=1:size(BaseEnds,1)
        %int = intersectEdges(BaseEnds(k,:), SegEnds(30:34,:));
        int = intersectEdges(BaseEnds(k,:), PatchEnds);
        crossind(k) = sum(sum(~isnan(int)))>0;
    end
else
    crossind = [];
end



[baselines, ib, ic] = unique(baselines,'rows');  %baselines(ic) = original baselines
tri_num = tri_num(ib);
BaseEnds = BaseEnds(ib,:);
if ~isempty(crossind)
    crossind = crossind(ib);
end


%unit directional vectors
Vec_unit = xy(baselines(:,1),:)-xy(baselines(:,2),:);
L = sqrt(Vec_unit(:,1).^2+Vec_unit(:,2).^2);
Vec_unit = Vec_unit./[L L];  



%dot velocity vectors into unit direction 
Vbase2 = Veast(baselines(:,2)).*Vec_unit(:,1) + Vnorth(baselines(:,2)).*Vec_unit(:,2);
Vbase1 = Veast(baselines(:,1)).*Vec_unit(:,1) + Vnorth(baselines(:,1)).*Vec_unit(:,2);

%propagate errors
Sig_base2 = sqrt(Vec_unit(:,1).^2.*Sigeast(baselines(:,2)).^2 + Vec_unit(:,2).^2.*Signorth(baselines(:,2)).^2);
Sig_base1 = sqrt(Vec_unit(:,1).^2.*Sigeast(baselines(:,1)).^2 + Vec_unit(:,2).^2.*Signorth(baselines(:,1)).^2);


Vbase = Vbase1 - Vbase2;
%propagate error
Sigbase = sqrt(Sig_base2.^2 + Sig_base1.^2);









