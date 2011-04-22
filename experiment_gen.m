function experiment_gen(exp_path, orig_hl_tpl_path, hl_basis, orig_tpl_path, n, p)
    mkdir(exp_path);
    tpl_name = 'tpl';
    tpl_path = [exp_path, '\\', tpl_name, '.tpl'];
    copyfile(orig_tpl_path, tpl_path);

    % FIXME: shouldn't need this additional template. This is workaround
    % for readfchk bug.
    hl_tpl_name = 'hl_tpl';
    hl_tpl_path = [exp_path, '\\', hl_tpl_name, '.tpl'];
    copyfile(orig_hl_tpl_path, hl_tpl_path);
    
    % Generate high-level model with fixed basis set

    frag_cfg = Fragment.defaultConfig();
    frag_cfg.template = hl_tpl_name;
    frag_cfg.basisSet = hl_basis;
    frag_cfg.par = [1.0 1.0 1.0 1.0 1.0];
    frag = Fragment(exp_path, frag_cfg);
    save([exp_path, '\\', 'frag.mat'], 'frag');

    % Generate low-level model with custom basis sets

    for i=1:p
        fraghat_cfg = Fragment.defaultConfig();
        fraghat_cfg.template = tpl_name;
        fraghat_cfg.basisSet = 'GEN';
        r = normrnd(1.0, 0.02, 1,3);
        fraghat_cfg.par = [r(1) r(2) r(2) r(3) r(3)];

        fraghats{i} = Fragment(exp_path, fraghat_cfg);
    end
    save([exp_path, '\\', 'fraghats.mat'], 'fraghats');

    % Generate environments and pair up adjacent ones

    env(1,2*n) = Environment;
    for i = 1:2*n
      env(1,i) = Environment.newCube([6,6,6],3);
    end
    env_pairs = cell(n, 1);
    for i = 1:n
      env_pairs{i} = {env(2*i - 1), env(2*i)};
    end

    save([exp_path, '\\', 'env_pairs.mat'], 'env_pairs');
end