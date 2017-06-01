function [w, wGrid, lambda, varPsi] = incValue(expDelta, expMu)

global n_grid n_period n_person n_ar;
global vars indx;

%vars.expIncValue: [n_period*n_product,n_person]
vars.num = (expDelta*indx.onesRowPerson) .* expMu;
vars.expIncValue = indx.sum * vars.num; 
w = log(vars.expIncValue);
vars.Wstack(n_ar+1:2*n_ar,1:n_person) = w(1:n_ar,1:n_person);
vars.Wflat = reshape(vars.Wstack, n_ar*n_person*2,1);
x = sparse(indx.iW, indx.jW, vars.Wflat);
y = reshape(w(2:n_period,1:n_person), n_ar*n_person,1);

lambda = (x'*x) \ (x'*y); % [2*n_person,1]
resid = reshape(y-x*lambda, n_ar, n_person);
varPsi = dot(resid,resid,1) / (n_ar-1);

max_w = max(w,[],1) + 2.575829*sqrt(varPsi);
min_w = min(w,[],1) - 2.575829*sqrt(varPsi);

%wGrid = repmat(min_w,n_grid,1) + (0:n_grid-1)'*(max_w-min_w)/(n_grid-1);
wGrid = bsxfun(@plus, min_w, (0:n_grid-1)'*(max_w-min_w)/(n_grid-1));
end

