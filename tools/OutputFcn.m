function varargout = OutputFcn(fn,varargin)
% OutputFcn works with SaveStreamData and Ode to record the solution at each
% time step to disk. This is much better than holding the solution in memory.
  [varargout{1:nargout}] = feval(fn,varargin{:});

function Init(saveFn,varargin)
  global ggof;
  ggof.saveEvery = process_options(varargin,'saveEvery',1);
  ggof.ctr = 0;
  ggof.ssd = SaveStreamData('Init',saveFn);
  
function stop = Fn(t,y,msg,varargin)
  global ggof;
  stop = 0;
  if(~isempty(msg))
    if(msg(1) == 'i')
      % Get the final time so we're sure to save the last point
      ggof.tf = t(end);
    else
      return;
    end
  end
  if(mod(ggof.ctr,ggof.saveEvery) == 0 || t == ggof.tf)
    ggof.ssd = SaveStreamData('Write',ggof.ssd,[t(1); y(:)]);
  end
  ggof.ctr = ggof.ctr + 1;
  
function Cleanup()
  global ggof;
  SaveStreamData('Finalize',ggof.ssd);
  clear global ggof;
  