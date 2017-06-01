function valfunOut = nfpBellman(valfunIn, wGrid, lambda, varPsi)

global tolBellman beta;

weight = transProb(wGrid, lambda, varPsi);
wGridFlat = wGrid(:);
dist = 1;
iter = 1;
while (dist > tolBellman && iter <=1500)
  valfunOut = log(exp(beta*weight*valfunIn - wGridFlat) + 1) + wGridFlat;
  diff = abs(valfunOut-valfunIn);
  dist = max(diff);
  valfunIn = valfunOut;
  iter = iter+1;
end

end

