function vec = mvp(xi,id,rs,cs,aslip_index)

x=zeros(size(aslip_index));
x(aslip_index) = xi;

vec = hm_mvp('mvp',id,x,rs,cs);

vec = vec(aslip_index);

end

