clear;

global iv dat mat par rnd vars indx ...
  n_product n_period n_person n_mc n_seed ...
  n_theta1 n_theta2 n_theta3 n_param ...
  n_inst n_grid n_ar n_obs ...
  n_choice dimValfun;

global beta tolBellman;

dat  = [];
mat  = [];
rnd  = [];
indx = [];
vars = [];
n_product = 50;
n_period  = 100;
n_person  = 50;
n_mc      = 20;
n_seed    = 5;

n_theta1 = 5;
n_theta2 = 3;
n_theta3 = 9;
n_inst   = 66;
n_param = 1+n_theta1+n_theta2+n_theta3;

n_grid = 50;
n_ar = n_period-1;
n_obs = n_period*n_product;
n_choice  = n_period;
dimValfun = n_person * n_grid;

tolBellman = 5e-15;

beta =    0.99;
beta0 =   6.00;
betaX =  [1.00, 1.00, 0.50];
alpha =   2.00;
gamma0 =  1.00;
gammaX = [0.20, 0.20, 0.10];
gammaZ =  1;
gammaW =  0.20;
gammaXi = 0.70;
rho0 =    0.10;
rhoZ =    0.95;
sdX =    [0.50, 0.50, 0.50];
sdXi =    1.00;
sdU =     0.01;
sdEta =   0.10;

rng(12345,'dsfmt19937');
filename = 'dat_dynblp4.txt';
prodIndx = 1:n_product;

par.theta1 = [beta0 betaX alpha]';
par.theta2 = [0.5; 0.5; 0.25];
boundUpper = 10*par.theta2;

cDataFile = fopen(filename, 'w');
cRandFile = fopen('rand.txt','w');
cSeedFile = fopen('seed.txt','w');
cValfFile = fopen('valfun.txt','w');

allocate();

for mcID = 1:n_mc,
  nu = randn(n_person,n_theta2)';
  x = zeros(n_product,n_period,n_theta1-2);
  for k=1:3,
    x(:,:,k) = sdX(k) * randn(n_product,n_period);
  end
  
  w   = randn(n_product,n_period);
  xi  = sdXi * randn(n_product,n_period);
  u   = sdU * randn(n_product,n_period);
  eta = sdEta * randn(n_product,n_period);
  
  z = zeros(n_product,n_period);
  z(:,1) = rho0 + rhoZ*8 + eta(:,1);
  for t=2:n_period,
    z(:,t) = rho0 + rhoZ*z(:,t-1) + eta(:,t);
  end
  p = gamma0 + gammaZ*z + gammaW*w + gammaXi*xi + u;
  for k=1:3,
    p = p + gammaX(k)*x(:,:,k);
  end
  % BLP inst
  iv.blp = zeros(n_product,n_period,size(x,3));
  for j=1:n_product,
    for k=1:size(x,3),
      iv.blp(j,:,k) = sum(x(prodIndx~=j,:,k));
      p(j,:) = p(j,:) - 0.1*iv.blp(j,:,k);
    end
  end
  
  iv.z = bsxfun(@plus, z(:), 0.5*randn(n_obs,3));
  iv.w = bsxfun(@plus, w(:), 0.5*randn(n_obs,3));
  xvec = reshape(x, n_obs,1,n_theta1-2);
  X = [ones(n_obs,1) squeeze(xvec) -p(:)];  
  A  = X(:,2:end-1);
  IV = [iv.z iv.w reshape(iv.blp,n_obs,size(x,3))];
  Z  = [X(:,1:end-1) A.^2 A.^3 IV IV.^2 IV.^3 prod(A,2) prod(IV,2) ...
    bsxfun(@times,A(:,1),IV) bsxfun(@times,A(:,2),IV) bsxfun(@times,A(:,3),IV)];
  
  delta = X*par.theta1 + xi(:);
  
  % Matrix constants
  %n_inst = size(Z,2);
  ZZ = Z' * Z;
  invPHI = ZZ \ eye(n_inst); %inv(mat.ZZ);
  XZ = X' * Z;
  XZPHIZXZ = (XZ*invPHI*(XZ)') \ (XZ*invPHI*Z');
  Xrc = X(:,[2,3,5]);
  
  valfunIn  = zeros(n_grid*n_person,1);
  
  mu = diag(par.theta2)*nu;
  expMu = exp(Xrc*mu);

  expDelta = exp(delta);
  [w, wGrid, lambda, varPsi] = incValue(expDelta, expMu);
  valfunOut = nfpBellman(valfunIn,wGrid,lambda,varPsi);
  EV = expectValfun(valfunOut,w,wGrid,lambda,varPsi);
  [ms, ms0] = choicePath(EV);
  
  % Store DATa to Matlab structs
  dat.ms(:,mcID) = ms;
  dat.ms0(:,mcID) = ms0;
  mat.X(:,:,mcID) = X;
  mat.Z(:,:,mcID) = Z;
  mat.invPHI(:,:,mcID) = invPHI;
  mat.XZPHIZXZ(:,:,mcID) = XZPHIZXZ;
  mat.Xrc(:,:,mcID) = Xrc;
  rnd.nu(:,:,mcID) = nu;
  for seed=1:n_seed
    rnd.seed(:,seed,mcID) = boundUpper .* rand(n_theta2,1);
  end
  
  % Export DATa to forMATted text files for C++
  cDataFile = fopen(filename, 'a');
  for t=1:n_period
    for j=1:n_product
      id = n_product*(t-1) + j;
      fprintf(cDataFile,'%5d %5d %5d %21.15f %21.15f %21.15f ',...
        mcID,t,j,ms(id), ms0(t),xi(j,t));
      for k=1:n_theta1
        fprintf(cDataFile,'%21.15f ',X(id,k));
      end
      for k=1:n_inst
        fprintf(cDataFile,'%21.15f ',Z(id,k));
      end
      fprintf(cDataFile,'\n');
    end
  end
  fclose(cDataFile);
  
  cRandFile = fopen('rand.txt','a');
  for n=1:n_person
    for k=1:n_theta2
      fprintf(cRandFile,'%21.15f ', nu(k,n));
    end
    fprintf(cRandFile,'\n');
  end
  fclose(cRandFile);
  
  cSeedFile = fopen('seed.txt','a');
  for seed=1:n_seed
    for k=1:n_theta2
      fprintf(cSeedFile,'%21.15f ',rnd.seed(k,seed,mcID));
    end
    fprintf(cSeedFile,'\n');
  end
  fclose(cSeedFile);
  
  cValfFile = fopen('valfun.txt','a');
  for n=1:n_person
    for g=1:n_grid
      fprintf(cValfFile,'%21.15f ', valfunOut(n_grid*(n-1)+g));
    end
    fprintf(cValfFile,'\n');
  end
  fclose(cValfFile);
  fprintf('MC %d complete\n', mcID);
  
end

save dynblp4.mat dat mat rnd indx boundUpper;
