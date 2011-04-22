function [test_risk] = experiment_run(exp_path, A, B)
    load([exp_path, '\\', 'frag.mat']);
    load([exp_path, '\\', 'fraghats.mat']);
    load([exp_path, '\\', 'env_pairs.mat']);
     
    p = size(fraghats, 1);
    n = size(env_pairs, 1);

    % Shuffle env pairs and split into train and test sets
    env_pairs = env_pairs(randperm(n));
    train_indices = arr(1:(n/2));
    test_indices = arr((n/2)+1:n);
    
    % Train and test
    aggfrag = AggregateFragment.train(frag, fraghats, train_indices);
    
    testlen = len(test_indices);
    test_risk = 0;
    for i=1:testlen
        index = test_indices{i};
        test_risk = test_risk + loss(frag, aggfrag, 2*index, 2*index + 1, A, B);
    end
end

