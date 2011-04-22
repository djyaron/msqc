function experiment_comp(exp_path)
    load([exp_path, '/', 'frag.mat']);
    load([exp_path, '/', 'fraghats.mat']);
    p = size(fraghats, 2);
    load([exp_path, '/', 'env_pairs.mat']);
    n = size(env_pairs, 1);

    % Perform all computations with Gaussian
    % this is easier due to the memoization
    for i=1:n
        env0 = env_pairs{i}{1};
        env1 = env_pairs{i}{2};
        
        frag.addEnv(env0);
        frag.addEnv(env1);
        
        for j=1:p
            fraghats{j}.addEnv(env0);
            fraghats{j}.addEnv(env1);
        end
    end
end

