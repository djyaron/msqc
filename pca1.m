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
hold off;