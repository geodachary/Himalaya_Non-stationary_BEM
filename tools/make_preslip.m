
load preslip
load pmpost_jess

%find patch on my fault that is closest to patch on Jess's fault
for k=1:size(pm,1)
    dist=sqrt((pm(k,6)-pmpost(:,6)).^2+(pm(k,7)-pmpost(:,7)).^2+(pm(k,3)-pmpost(:,3)).^2);
    [dummy index]=min(dist);
    preslip_new(k)=preslip(index);
 end

% index = pm(:,7)<-10 & pm(:,3)>8 & pm(:,3)<15;
% preslip_new(index)=10^-5;