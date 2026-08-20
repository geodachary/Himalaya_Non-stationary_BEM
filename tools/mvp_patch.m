function vec = mvp(xi,id,rs,cs,aslip_index,fixpatch_source,fixpatch_rec)

x=zeros(size(aslip_index));
x(aslip_index) = xi;

vec = hm_mvp('mvp',id,x,rs,cs);


%fix problem source by setting problem receiver to zero due to source
if ismember(fixpatch_rec,rs)
    if ismember(fixpatch_source,cs)

    cs(cs==fixpatch_source) = [];
    vec_patch = hm_mvp('mvp',id,x,fixpatch_rec,cs);
 
    vec(fixpatch_rec) = vec_patch(fixpatch_rec);
         
    end
end

vec = vec(aslip_index);

end

