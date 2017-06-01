function weight = transProb(wGrid, lambda, varPsi)

global n_grid n_person;
global indx;

lambdaMat = reshape(lambda, 2,n_person);
% lambda0   = repmat(lambdaMat(1,:), n_grid,1);
% lambda1   = repmat(lambdaMat(2,:), n_grid,1);
% wLambda   = lambda0 + lambda1.*wGrid;
lambdaW     = bsxfun(@times,lambdaMat(2,:),wGrid);
condMean    = bsxfun(@plus,lambdaMat(1,:),lambdaW);
condMeanExp = kron(condMean, ones(1,n_grid));
% wGridMat  = repmat(reshape(wGrid, 1, n_grid*n_person), n_grid,1);
% wMatExp   = wGridMat - condMeanExp;
wMatExp = bsxfun(@minus, reshape(wGrid,1,n_grid*n_person), condMeanExp);

% varPsiExp = repmat(kron(varPsi, ones(1,n_grid)), n_grid,1);
% numWeight = exp(-0.5*wMatExp .* wMatExp ./ varPsiExp);
varPsiExp = kron(varPsi, ones(1,n_grid));
numWeight = exp(-0.5*wMatExp .* bsxfun(@rdivide,wMatExp,varPsiExp));
sumWeight = numWeight * indx.sumW;
denWeight = kron(1./sumWeight, ones(1,n_grid));
weightMat = numWeight .* denWeight;
weight    = sparse(indx.iWeight,indx.jWeight,weightMat(:));

end

