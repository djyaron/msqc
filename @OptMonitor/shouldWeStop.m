function stop = shouldWeStop(obj)

values = obj.values{end};
niter = values.iteration;
if (niter > obj.maxIter)
   stop = 1;
   obj.stoppingCriterion = 'exceeded max iterations';
   return
end
% don't want to terminate before 4 successful steps
needed = 5;
ngood = length(obj.etest);
if (ngood == 0)
   disp(['train ',num2str(obj.tokcal(obj.values{end}.residual))]);
   stop = 0;
   return;
elseif (ngood < needed)
   disp(['train ',num2str(obj.tokcal(obj.values{end}.residual)), ...
      ' test ',num2str(obj.tokcal(obj.etest{end}))]);
   stop = 0;
   return;
end
stop = 0;
if (obj.ediff > 0)
   % create list of changes in energy from past iterations
   improvement = zeros(needed-1,1);
   rcurr = obj.tokcal(obj.etest{end});
   for i=1:(needed-1)
      improvement(i) = obj.tokcal(obj.etest{end-i})-rcurr;
   end
   disp(['train ',num2str(obj.tokcal(obj.values{end}.residual)), ...
      ' test ', num2str(rcurr), ...
      ' improvements ',num2str(improvement(:)')]);
   if (max(improvement) < obj.ediff)
      stop = 1;
      obj.stoppingCriterion = ['Ediff converged ',... 
         num2str(improvement(:)')];
   end
end
% if (obj.xdiff > 0)
%    % create list of changes in energy for past 4 iterations
%    n=4;
%    ed = zeros(4,1);
%    for i=1:n
%       ed(i) = norm(obj.param{end} - obj.param{end-i});
%    end
%    disp(['checking xdiff ',num2str(ed(:)')]);
%    check = max(abs(ed));
%    if (check < obj.xdiff)
%       stop = 1;
%       obj.stoppingCriterion = ['Xdiff converged ',num2str(ed(:)')];
%    end
% end   
end

