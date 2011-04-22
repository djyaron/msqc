function [test_risk] = experiment_run(exp_path, A, B)
    load([exp_path, '\\', 'frag.mat']);
    frag.loadAllEnv;
    load([exp_path, '\\', 'fraghats.mat']);
    p = size(fraghats, 2);
    for i=1:p
      fraghats{i}.loadAllEnv;
    end
    load([exp_path, '\\', 'env_pairs.mat']);
    n = size(env_pairs, 1);

    % Shuffle env pairs and split into train and test sets
    shuffle = randperm(n);
    train_indices = shuffle(1:(n/2));
    test_indices = shuffle((n/2)+1:n);
    
    % Train and test
    aggfrag = AggregateFragment.train(frag, fraghats, train_indices, A, B);
    
    testlen = size(test_indices, 2);
    test_risk = 0;
    for i=1:testlen
        i0 = 2 * test_indices(i) - 1;
        i1 = 2 * test_indices(i);
        test_risk = test_risk + loss(frag, aggfrag, i0, i1, A, B);
    end
    test_risk = test_risk / testlen;
end

