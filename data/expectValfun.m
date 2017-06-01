function EV = expectValfun(valfunIn, w, wGrid, lambda, varPsi)

global n_person n_grid;
global indx;

lambdaMat   = reshape(lambda, 2,n_person);
% lambda0   = repmat(lambdaMat(1,:), n_period,1);
% lambda1   = repmat(lambdaMat(2,:), n_period,1);
% condMean  = lambda0 + lambda1.*w; % [n_period,n_person]
lambdaW     = bsxfun(@times,lambdaMat(2,:),w);
condMean    = bsxfun(@plus,lambdaMat(1,:),lambdaW);
condMeanExp = kron(condMean, ones(1,n_grid));
% wGridMat    = repmat(reshape(wGrid, 1, n_grid*n_person), n_period,1);
% wMatExp     = wGridMat - condMeanExp;
wMatExp     = bsxfun(@minus, reshape(wGrid, 1, n_grid*n_person), condMeanExp);

% varPsiExp = repmat(kron(varPsi, ones(1,n_grid)), n_period,1);
% numWeight = exp(-0.5*wMatExp .* wMatExp ./ varPsiExp);
varPsiExp = kron(varPsi, ones(1,n_grid));
numWeight = exp(-0.5*wMatExp .* bsxfun(@rdivide,wMatExp,varPsiExp));
sumWeight = numWeight * indx.sumW;
denWeight = kron(1./sumWeight, ones(1,n_grid));
weightMat = numWeight .* denWeight;
weight    = sparse(indx.iWeight2,indx.jWeight2,weightMat(:));

% lambdaMat = reshape(lambda, 2,n_person);
% lambda0   = repmat(lambdaMat(1,:), n_period,1);
% lambda1   = repmat(lambdaMat(2,:), n_period,1);
% 
% condMean    = lambda0 + lambda1.*w; % [n_period,n_person]
% condMeanExp = kron(condMean, ones(1,n_grid));
% wGridMat    = repmat(reshape(wGrid, 1, n_grid*n_person), n_period,1);
% wMatExp     = wGridMat - condMeanExp;
% 
% varPsiExp = repmat(kron(varPsi, ones(1,n_grid)), n_period,1);
% numWeight = exp(-0.5*wMatExp .* wMatExp ./ varPsiExp);
% sumWeight = numWeight * indx.sumW;
% denWeight = kron(1./sumWeight, ones(1,n_grid));
% weightMat = numWeight .* denWeight;
% weight    = sparse(indx.iWeight2,indx.jWeight2,weightMat(:));

% wMatExp = repmat(wGrid',n_period,1) - repmat(condMean,1,n_grid);
% numWeight = exp(-0.5 * wMatExp .* wMatExp / varPsi);
% sumWeight = sum(numWeight,2);
% denomWeight = 1./sumWeight; %bsxfun(@rdivide,1,sumWeight);
% weightMat = numWeight .* repmat(denomWeight,1,n_grid);
EV = weight * valfunIn;

end

