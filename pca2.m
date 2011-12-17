%%
clear classes;
load(['data\ethane\ethane4\data0.mat']);
%%
clear classes;
load('data\ethane\ethane3\data.mat');
%% set up the pca data
dataSize = size( HL );
figure( 67 );
hold on;
for ifrag = 1:dataSize( 1 )
    frag1 = HL{ifrag,1};
    nenv = frag1.nenv;
    
    % Find the number of basis functions for each atom in the fragment.
    nbasis = zeros(1,frag1.natom);
    for iatom = 1:frag1.natom
        ind = find(frag1.basisAtom == iatom);
        nbasis(iatom) = size(ind,1);
    end
    nbasisSq = nbasis .^ 2;
    
    % Set up the PCA input matrix for each atom in each environment.
    % Density data for each environment is reduced to a single row.
    % Data for each atom is stored sequentially.
    X = zeros(nenv, sum(nbasisSq));
    for iatom = 1:frag1.natom
        % Basis functions for the given atom
        ind = find(frag1.basisAtom == iatom);
        
        % Find the indices within the row for that atom
        if iatom == 1
            x1 = 1;
        else
            x1 = sum( nbasisSq(1,1:iatom-1) ) + 1;                
        end
        x2 = sum( nbasisSq(1,1:iatom) );
        
        for ienv = 1:nenv
            % Add the density data to the input matrix.
            den = frag1.density(ienv);
            den1 = den(ind,ind);
            X(ienv,x1:x2) = reshape(den1, [1,nbasis(iatom)^2]);
        end
    end
    
    % Subtract off the means (does PCA automatically do this?)
    means = mean(X,1);
    for ienv=1:nenv
        X(ienv,:) = X(ienv,:)-means;
    end
    
    % Run PCA!. 'econ' removes zeros from the 'latent' result.
    [coef, score, latent] = princomp(X, 'econ');
    
    % Plot the covariance data on a logarithmic scale.
    LnLatent = log(latent);
    %cutoff = find( latent < 10^-25, 1, 'first' );
    %if size( cutoff, 1 ) > 0
    %    LnLatent = LnLatent( 1:cutoff, : );
    %end
    plot( LnLatent, 'b' );
end
hold off;