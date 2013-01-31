classdef OptMonitor < handle
   properties
      ftrain    % fitme for training, need pointer to update context
      ftest     % fitme object for monitoring convergence
      param     % cell array of params
      values    % cell array of optimValues
      state     % cell array of states ('init','iteration' etc)
      accepted  % list of accepted iterations
      etest     % cell array of test error on accepted iterations
      lsqOutput % all outputs returned by lsqnonlin
      maxIter % 
      ediff   %
      stoppingCriterion % string holding information on termination
      plotNum % if nonzero, plots to this figure number
      updateContext % update context after each successful iteration
   end
   methods (Static)
      function res = tokcal(errs)
         res = norm(errs)*627.509/sqrt(length(errs));
      end
   end
   methods
      function res = OptMonitor(ftrain,ftest,maxIter, ediff)
         if (nargin < 3)
            maxIter = 50;
         end
         if (nargin < 4)
            ediff = 0.01;
         end
         res.ftrain = ftrain;
         res.ftest = ftest;
         res.param = cell(0,0);
         res.values =  cell(0,0);
         res.state = cell(0,0);
         res.etest = cell(0,0);
         res.maxIter = maxIter;
         res.ediff   = ediff;
         res.accepted = [];
         res.plotNum = 0;
         res.updateContext = 0;
      end
      function stop = toCall(obj,x,optimValues,state,varargin)
         %disp('inside OptMontor.toCall');
         if (strcmpi(state,'interrupt'))
            % nothing has change in x or optimValues, but we are being
            % given the chance to terminate if we like
            stop = 0;
         elseif (strcmpi(state,'init'))
            figure(obj.plotNum);
            clf;
            stop = 0;
         else
            obj.param{end+1}  = x;
            obj.values{end+1} = optimValues;
            obj.state{end+1}  = state;
            if (length(varargin) > 1)
               disp('warning, extra info being passed to OptMonitor.toCall');
            end
            if (obj.plotNum)
               figure(obj.plotNum);
               hold on;
               plot(optimValues.iteration, ...
                  obj.tokcal(optimValues.residual),'bo');
            end
            % Did the iteration succeed
            if (length(obj.values) > 1)
               if (abs(obj.values{end}.resnorm - ...
                     obj.values{end-1}.resnorm) > 1e-12)
                  obj.accepted(end+1) = length(obj.values);
                  obj.etest{end+1} = obj.ftest.err(x);
                  if (obj.plotNum)
                     figure(obj.plotNum);
                     hold on;
                     plot(optimValues.iteration, ...
                        obj.tokcal(obj.etest{end}),'ro');
                  end
               end
               if (obj.updateContext)
                  obj.ftrain.updateContext;
                  obj.ftest.updateContext;
               end
            end
            stop = obj.shouldWeStop;
         end
      end
      function storeoutput(obj, pt,resnorm,residual, ...
            exitflag,output,lambda,jacobian)
         obj.lsqOutput.pt = pt;
         obj.lsqOutput.resnorm  = resnorm;
         obj.lsqOutput.residual = residual;
         obj.lsqOutput.exitflag = exitflag;
         obj.lsqOutput.output   = output;
         obj.lsqOutput.lambda   = lambda;
         obj.lsqOutput.jacobian = jacobian;
      end
      function [y,x] = trainError(obj)
         niter = length(obj.values);
         x = zeros(niter,1);
         y = zeros(niter,1);
         for i = 1:niter
            optimValues = obj.values{i};
            x(i) = optimValues.iteration;
            y(i) = obj.tokcal(optimValues.residual);
         end
      end
      function [y,x] = testError(obj)
         niter = length(obj.etest);
         x = zeros(niter,1);
         y = zeros(niter,1);
         for i = 1:niter
            igood = obj.accepted(i);
            optimValues = obj.values{igood};
            x(i) = optimValues.iteration;
            y(i) = obj.tokcal(obj.etest{i});
         end     
      end
      function [y,x] = getPardiff(obj)
         niter = length(obj.values);
         x = zeros(niter-1,1);
         y = zeros(niter-1,1);
         for i = 2:niter
            optimValues = obj.values{i};
            x(i-1) = optimValues.iteration;
            y(i-1) = norm(obj.param{i}-obj.param{i-1});
         end
      end
      function plotErrors(obj,figNum)
         [y,x] = obj.trainError;
         figure(figNum)
         plot(x,y,'bo');
         [y,x] = obj.testError;
         hold on;
         plot(x,y,'ro');
      end
   end
end

