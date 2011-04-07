classdef AggregateFragment < handle
    properties
        weights
        frags
    end
    
    methods(Static)
        function af = train(frag_cfg, fraghat_cfgs, tplpath, experiments)
            fraghats = cellfun(@(c) {tmpfrag(c, tplpath)}, fraghat_cfgs);
            m = size(fraghats, 2);
            % FIXME: memoization
            w = repmat(1/m, m); 
            af = AggregateFragment(w, fraghats);
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

