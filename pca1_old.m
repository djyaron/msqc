<<<<<<< HEAD
%%
clear classes;
%root = 'c:\dave\apoly\msqc\';
load(['data\ethane\ethane2\data3.mat']);
%%
clear classes;
load('data\ethane\ethane3\data.mat');
%% set up the pca data
figure( 67 );
hold on;
dataSize = size( HL );
for ifrag = 1:dataSize( 1 )
    frag1 = HL{ifrag,1};
    nenv = frag1.nenv;
    ind = find(frag1.basisAtom == 1);
    nbasis = size(ind,1);
    X = zeros(nenv, nbasis^2);
    for ienv = 1:nenv
      den = frag1.density(ienv);
      den1 = den(ind,ind);
      X(ienv,:) = reshape(den1, [1,nbasis^2]);
    end
    means = mean(X,1);
    for ienv=1:nenv
        X(ienv,:) = X(ienv,:)-means;
    end
    [coef, score, latent] = princomp(X, 'econ');
    LnLatent = log(latent);
    cutoff = find( latent < 10^-25, 1, 'first' );
    if size( cutoff, 1 ) > 0
        LnLatent = LnLatent( 1:cutoff, : );
    end
    plot( LnLatent, 'b' );
end
=======
clear classes;
%root = 'c:\dave\apoly\msqc\';
load(['data\ethane\ethane2\data3.mat']);
%% set up the pca data
dataSize = size( HL );
nenv = 100;
ic = 0;
X = zeros(7*nenv, nbasis^2);
for ifrag = 1:7
    frag1 = HL{ifrag,1};
    nenv =frag1.nenv;
    ind = find(frag1.basisAtom == 1);
    nbasis = frag1.nbasis;
    for ienv = 1:nenv
        den = frag1.density(ienv);
        %       den1 = den(ind,ind);
        ic = ic+1;
        X(ic,:) = reshape(den, [1,nbasis^2]);
    end
 
end

   Y = reshape(corrEnv', [1,700]);
    %     means = mean(X,1);
%     for ienv=1:nenv
%         X(ienv,:) = X(ienv,:)-means;
%     end
    [coef, score, latent] = princomp(X);
    figure( 67 );
    hold on;
    plot( latent );
    
    %solve score*beta = Y
    scorep = score(:,1:20);
    scorep = [ones(700,1) scorep];
    beta =scorep\Y';
    Ypred = scorep*beta;
    
    figure(68)
    plot(Ypred,Y,'r.')
>>>>>>> c5f65676a3b271e1d70f4b8e5152d64e978167f3
hold off;