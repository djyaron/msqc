classdef AggregateFragment < handle
    properties
        weights
        frags
    end
    
    methods(Static)
        function af = train(frag, fraghats, ep_indices, A, B)
            p = size(fraghats, 2);
            rate = 0.1;
            T = len(ep_indices);
            w = repmat(1/p, T, p); 
            x = zeros(p);
            
            for t=1:T
                i0 = ep_indices{t}{1};
                i1 = ep_indices{t}{2};
                
                y = Delta(frag, i0, i1, A, B);
                for j=1:p
                  x(j) = Delta(fraghats{j}, i0, i1, A, B);
                end
                yhat = w(t) . x;
                % Absolute loss
                sgn = sign(y - yhat);
                g = rate * sgn * x; 
                
                % Gradient descent
                % w(t+1) = w(t) - g;
                % Exponentiated gradient descent (Warmuth 97)
                r = exp(g);
                w(t+1) = r .* w(t) / (r.w(t));  
            end
            
            af = AggregateFragment(sum(w), fraghats);
        end
    end
    
    methods 
        function c = comb(obj, l)
          c = dot(obj.weights, l);
        end
    end
    
    methods
        function af = AggregateFragment(w, f)
          af.weights = w;
          af.frags = f;
        end
        
        function h = H1(obj)
            h = comb(cellfun(@(f) f.H1, obj.frags));
        end
        
        function h = H1Env(obj)
            h = comb(cellfun(@(f) f.H1Env, obj.frags));
        end
        
        function h = H2(obj)
            h = comb(cellfun(@(f) f.H2, obj.frags));
        end
        
        function d = density(obj, env)
            d = comb(cellfun(@(f) f.density(env), obj.frags));
        end
        
        function density2p(obj, env)
            comb(cellfun(@(f) f.density2p(env), obj.frags));
        end
    end
end

