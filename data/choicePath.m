function [ms,ms0] = choicePath(EV)

global n_period n_product n_person;
global beta;
global vars;

expValue0 = reshape(exp(beta*EV), n_period, n_person);
denom = 1./(vars.expIncValue + expValue0);
denExpand = kron(denom, ones(n_product,1));
choiceProb = vars.num .* denExpand; % [n_period*n_product,n_person]

nochoiceProb = expValue0.*denom;
surviveProb = [ones(1,n_person); cumprod(nochoiceProb(1:n_period-1,1:n_person))];
surviveProbExpand = kron(surviveProb, ones(n_product,1));

ms = (choiceProb .* surviveProbExpand) * ones(n_person,1) / n_person;
ms0 = (nochoiceProb .* surviveProb) * ones(n_person,1) / n_person;

end

