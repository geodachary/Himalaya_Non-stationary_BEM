function varargout = CalcGds_tri_O(fn,varargin)
% Rearranged code from make_stressG_ss.m.
  [varargout{1:nargout}] = feval(fn,varargin{:});

function Make(el,nd,rerr,savefn)
  g = CalcGds_tri('Init',el,nd);
  CalcGds_tri_O('Compress',g,rerr,savefn);
  
function g = Init(el,nd)
    
    patch_stuff=make_triangular_patch_stuff(el,nd);


    centroids=patch_stuff.centroids_faces;
    normal=patch_stuff.normal_faces;
    strikevec=patch_stuff.strikevec_faces;
    dipvec=patch_stuff.dipvec_faces;

   


  g.centroids=centroids;
  g.strikevec=strikevec;
  g.normal=normal;
  
  g.mu = 3e10;
  g.dipvec = dipvec;

  g.nd = nd;
  g.el = el;
  
function Gds = Fn(g,rs,cs)



  Gds = zeros(length(rs),length(cs));
  for(i = 1:length(cs))
    k = cs(i);
    
     temp1{1}=g.nd;
     temp2{1}=g.el(k,:);    

     ss = 0;
     ds = -1; %negative so that self-stress is negative (stress drop)

      Strain = CalcTriStrains(g.centroids(rs,1), g.centroids(rs,2), g.centroids(rs,3), temp1{1}(temp2{1},1), temp1{1}(temp2{1},2), temp1{1}(temp2{1},3), 0.25, ss, 0, ds);    
       S = StrainToStress(Strain, g.mu, g.mu);
   

     
     normal = g.normal(rs,:);
%     %calculate tractions


    T1=S.xx.*normal(:,1)+S.xy.*normal(:,2)+S.xz.*normal(:,3);
    T2=S.xy.*normal(:,1)+S.yy.*normal(:,2)+S.yz.*normal(:,3);
    T3=S.xz.*normal(:,1)+S.yz.*normal(:,2)+S.zz.*normal(:,3);

    
    %component along dipvec
     dipvec = g.dipvec(rs,:);
     Sd_ds=T1.*dipvec(:,1)+T2.*dipvec(:,2)+T3.*dipvec(:,3);

     


     Gds(:,i) = Sd_ds;
       
   
     
  %  end
  
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
  end
  
  
function bs = Compress(g,rerr,savefn)


  G = @(rs,cs)Fn(g,rs,cs);
  
   %X is 3xN point cloud
   X = (g.centroids)';
   
 
    [bs p q] = hmcc_hd(X,X);

    t = tic();
    n = length(p);
    eBfro = hm('EstBfroLowBndCols',G,p,q,1:n,1:n);
    savefn_tmp = [savefn '.tmp'];
    hm('Compress',bs,p,q,G,1e-5*rerr,savefn_tmp,'Bfro',eBfro);
    et = toc(t);
    [id hnnz] = hm_mvp('init',savefn_tmp);
    in = hm('HmatInfo',savefn_tmp);
    cf = in.m*in.n/hnnz;
    Bfro = sqrt(hm_mvp('fronorm2',id));
    [ce re] = hm('CondEstFro',id,struct('reltol',0.05));
    hm_mvp('cleanup',id);
    fprintf(1,'et: %1.1f  cf: %1.1f  ||B||_F: %1.1e %1.1e  ce: %1.1e\n',...
	    et,cf,eBfro,Bfro,ce);
    
    if (false)
      t = tic();
      hm('Compress',bs,p,q,G,rerr,savefn,...
	 'old_hmat_fn',savefn_tmp,'tol_method','BinvFro','BinvFro',ce/Bfro);
      [id hnnz] = hm_mvp('init',savefn);
      hm_mvp('cleanup',id);
      cf = in.m*in.n/hnnz;
      fprintf(1,'et: %1.1f  cf: %1.1f\n',toc(t),cf);
    else
      system(sprintf('mv %s %s',savefn_tmp,savefn));
    end
 
  
function v = uvec(v)
  v = v./repmat(sqrt(sum(v.^2)),size(v,1),1);
