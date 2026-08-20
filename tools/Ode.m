function varargout = Ode(fn,varargin)
% Explicit 2-3 pair (same as Matlab's). ode('Pair23') is meant to replace
% Matlab's ode23 function when a large problem is solved. It does not save more
% in memory than necessary. A user-supplied function is given the solution at
% each time step to save by the method desired.
%
% o = ode('Options',...): Same as 'odeset'. Use
%     stop = OutputFcn(t,y,msg,...)
% to record the solution.
%
% ode('Pair23',odeFcn,tspan,y0,opts,...): Same call sequence as ode23. However,
% there is one more output. The second output is set to 0 nominally. To signal
% an internal error in odeFcn and thereby request a shorter time step, set the
% second output to 1.

  [varargout{1:nargout}] = feval(fn,varargin{:});
  
function o = Options(varargin)
% Can use this or just call odeset.
  [o.RelTol o.AbsTol o.OutputFcn o.InitialStep] =...
      process_options(varargin,'RelTol',1e-3,'AbsTol',1e-6,...
		      'OutputFcn',[],'InitialStep',1e-3);
  
function Pair23(odeFcn,tspan,y0,opts,varargin)
% Clone of Matlab's ode23. Cite
%   P. Bogacki, L.F. Shampine, "A 3(2) Pair of Runge-Kutta Formulas",
%   Appl. Math Lett. 2(4), 321-325, 1989.

  threshold = opts.AbsTol(:)/opts.RelTol;
  tdir = sign(diff(tspan));
  n = length(y0);
  
  hof = ~isempty(opts.OutputFcn);
  if(hof) of = opts.OutputFcn; end

  f = zeros(n,4);
  t = tspan(1);
  absh = opts.InitialStep;

  
  [f(:,1) serr] = feval(odeFcn,t,y0,varargin{:});
  y = y0;
  
  if(hof)
    feval(of,tspan,y,'init',varargin{:});
  end
  
  % Taken from ode23.m
  pow = 1/3;
  A = [1/2, 3/4, 1];
  B = [1/2         0               2/9
       0           3/4             1/3
       0           0               4/9
       0           0               0  ]; 
  E = [-5/72; 1/12; 1/9; -1/8];
  
  while(t < tspan(end)) % Loop to integrate
    hmin = 16*eps(t);
    nofailed = true;
    while(1) % Loop for one step
      h = tdir*absh;
      if(t + h > tspan(end))
	h = tspan(end) - t;
	absh = abs(h);
      end
    
      hB = h*B;
      [f(:,2) serr] = feval(odeFcn,t+h*A(1),y+f*hB(:,1),varargin{:});
      if(~serr)
	[f(:,3) serr] = feval(odeFcn,t+h*A(2),y+f*hB(:,2),varargin{:});
	if(~serr)
	  yt = y + f*hB(:,3);
	  [f(:,4) serr] = feval(odeFcn,t+h*A(3),yt,varargin{:});
	  if(~serr)
	    % h*f*E is the difference between the 3rd- and 2nd-order methods
	    err = absh*norm((f*E)./max(max(abs(y),abs(yt)),threshold),inf);
	  end
	end
      end
      if(serr)
	err = inf;
	nofailed = false;
      end
	    
      if(err > opts.RelTol)
	if(absh <= hmin)
	  if(hof)
	    feval(of,t,y,'tolfail',varargin{:});
	  end
	  warning(sprintf('ode: Failure at t=%e. Unable to meet tol.\n',t));
	  return;
	end
	
	if(nofailed)
	  nofailed = false;
	  absh = max(hmin,absh*max(0.5,0.8*(opts.RelTol/err)^pow));
	else
	  absh = max(hmin,0.5*absh);
	end
      else
	break;
      end
    end

    if(nofailed)
      fac = 0.8*(opts.RelTol/err)^pow;
      absh = min(5,fac)*absh;
    end
    
    t = t + h;
    y = yt;
    f(:,1) = f(:,4);
    if(hof)
      stop = feval(of,t,y,'',varargin{:});
      if(stop) return; end
    end
  end
