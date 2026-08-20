function varargout = CalcG_rake(fn,varargin)
% Rearranged code from make_stressG_ss.m.
  [varargout{1:nargout}] = feval(fn,varargin{:});
 
function Make(el,nd,rake,rerr,savefn)
  g = CalcG_rake('Init',el,nd,rake);
  CalcG_rake('Compress',g,rerr,savefn);
  
function g = Init(el,nd,rake)
    
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

  
  %create rake vector as a rotation of strike vector about normal
    
    R = zeros(3,3);
    for k=1:size(g.normal,1)
        n = g.normal(k,:);  %negative so normal points up
        a = rake(k)*pi/180;  
        R(1,:) = [cos(a)+n(1)^2*(1-cos(a))          n(1)*n(2)*(1-cos(a))-n(3)*sin(a)  n(1)*n(3)*(1-cos(a))+n(2)*sin(a)];
        R(2,:) = [n(1)*n(2)*(1-cos(a))+n(3)*sin(a)  cos(a)+n(2)^2*(1-cos(a))          n(2)*n(3)*(1-cos(a))-n(1)*sin(a)];
        R(3,:) = [n(1)*n(3)*(1-cos(a))-n(2)*sin(a)  n(2)*n(3)*(1-cos(a))+n(1)*sin(a)  cos(a)+n(3)^2*(1-cos(a))];

        rakevec(k,:) = (R*g.strikevec(k,:)')';
    end
  
  g.rakevec = rakevec;

  
  
  g.nd = nd;
  g.el = el;
  g.rake = rake;
  
  
function Grake = Fn(g,rs,cs)



  Grake = zeros(length(rs),length(cs));
  for(i = 1:length(cs))
    k = cs(i);
    
     nd=g.nd;
     tri=g.el(k,:);   
     rake = g.rake(k);
     ss = cos(rake*pi/180);
     ds = sin(rake*pi/180);
     
    %[U, D, S] = tridisloc3d(g.centroids(rs,:)', temp1{1}', temp2{1}', [ss dd 0]', 1, .25);
     
     %note: negative signs have been carefully checked for consistency with definitation of rake
    [d,Strain] = CalcTriStrains_O(g.centroids(rs,1), g.centroids(rs,2), g.centroids(rs,3), nd(tri,1), nd(tri,2), nd(tri,3), 0.25, -ss, 0, -ds);   
    S = StrainToStress(Strain, 1, 1);


    %%triangular dislocation solution of Nikkhoo (below) is slower than
    %%Owen's code (above)
    
      % Calculated stress tensor components in EFCS. The six columns of Stress 
    % are Sxx, Syy, Szz, Sxy, Sxz and Syz, respectively. The stress components 
    % have the same unit as Lame constants.
        
     %note: signs on ss and ds have been carefully checked for consistency with definitation of rake
     %[S,Strain]=TDstressHS(g.centroids(rs,1), g.centroids(rs,2), g.centroids(rs,3),  nd(tri(1),:),  nd(tri(2),:),  nd(tri(3),:),-ss,ds,0,1,1);

  
     normal = g.normal(rs,:);
%     %calculate tractions

%if using Nikkhooo
    %T1=S(:,1).*normal(:,1)+S(:,4).*normal(:,2)+S(:,5).*normal(:,3);
    %T2=S(:,4).*normal(:,1)+S(:,2).*normal(:,2)+S(:,6).*normal(:,3);
    %T3=S(:,5).*normal(:,1)+S(:,6).*normal(:,2)+S(:,3).*normal(:,3);

%if using Owen
     T1=S.xx'.*normal(:,1)+S.xy'.*normal(:,2)+S.xz'.*normal(:,3);
     T2=S.xy'.*normal(:,1)+S.yy'.*normal(:,2)+S.yz'.*normal(:,3);
     T3=S.xz'.*normal(:,1)+S.yz'.*normal(:,2)+S.zz'.*normal(:,3);

    
    %component of traction in rake direction
    rakevec = g.rakevec(rs,:);
    Srake = T1.*rakevec(:,1) + T2.*rakevec(:,2) + T3.*rakevec(:,3);
   
     


     Grake(:,i) = -Srake; %negative sign so that stress drops with slip
  
       
   
     
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
