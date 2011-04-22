classdef AggregateFragment < handle
    properties
        weights
        frags
        
        % XXX need uniform interface to Fragment and AggregateFragment
        H1
        H1Env
        H2
    end
    
    methods(Static)
        function af = train(frag, fraghats, ep_indices, A, B)
            p = size(fraghats, 2);
            rate = 0.1;
            T = size(ep_indices, 2);
            w = repmat(1/p, T, p); 
            x = zeros(1,p);
            
            for t=1:T
                i0 = 2 * ep_indices(t) - 1;
                i1 = 2 * ep_indices(t);
                
                y = Delta(frag, i0, i1, A, B);
                for j=1:p
                  x(j) = Delta(fraghats{j}, i0, i1, A, B);
                end
                yhat = dot(w(t,:), x);
                % Absolute loss
                sgn = sign(y - yhat);
                g = rate * sgn * x; 
                
                % Gradient descent
                % w(t+1,:,) = w(t,:) - g;
                % Exponentiated gradient descent (Warmuth 97)
                r = exp(g);
                w(t+1,:) = r .* w(t) / dot(r, w(t,:)');  
            end
            
            wbar = sum(w(2:T+1,:)) / T;
            af = AggregateFragment(wbar, fraghats);
        end
    end
    
    methods 
        function c = comb(obj, l)
          c=0;
          for i=1:size(obj.weights)
            c = c + obj.weights(i) * l{i};
          end
        end
    end
    
    methods
        function af = AggregateFragment(w, f)
          af.weights = w;
          af.frags = f;
          
          af.H1 = af.comb(cellfun(@(f) f.H1, af.frags, 'UniformOutput', false));
          af.H1Env = af.comb(cellfun(@(f) f.H1Env, af.frags, 'UniformOutput', false));
          af.H2 = af.comb(cellfun(@(f) f.H2, af.frags, 'UniformOutput', false));
        end
             
        function d = density(obj, env)
            d = obj.comb(cellfun(@(f) f.density(env), obj.frags, 'UniformOutput', false));
        end
        
        function d = density2p(obj, env)
            d = obj.comb(cellfun(@(f) f.density2p(env), obj.frags, 'UniformOutput', false));
        end
    end
end

