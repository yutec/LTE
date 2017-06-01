function [] = allocate()

global n_period n_product n_person n_grid n_ar n_obs  ...
  n_theta1 n_theta2 n_inst n_mc n_seed;
global DAT MAT RND indx vars;

% Variables
vars.num         = zeros(n_obs,n_person);
vars.expIncValue = zeros(n_period,n_person);
vars.Wstack      = ones(n_ar*2,n_person);
vars.Wflat       = ones(n_ar*2*n_person,1);

% Index objects for vectorized operation
indx.onesColProduct = ones(n_product,1);
indx.expandRowTbyJ  = kron((1:n_period)',ones(n_product,1));
indx.onesRowPerson  = ones(1,n_person);
indx.onesColPerson  = indx.onesRowPerson';
indx.onesRowGrid    = ones(1,n_grid);
indx.onesColGrid    = indx.onesRowGrid';
%indx.oo = ones(1,n_person); % expand mean value
indx.sum = kron(speye(n_period), ones(1,n_product)); %for incValue
indx.sumW = kron(speye(n_person), ones(n_grid,1));
indx.iW = reshape((1:n_ar*n_person), n_ar,n_person);
indx.iW = reshape(kron(indx.iW,ones(1,2)), n_ar*n_person*2,1);
indx.jW = reshape(kron(1:n_person*2,ones(n_ar,1)), n_ar*n_person*2,1);
indx.iWeight = reshape(1:n_grid*n_person, n_grid,n_person);
indx.iWeight = kron(indx.iWeight, ones(1,n_grid));
indx.iWeight = indx.iWeight(:);
indx.jWeight = kron((1:n_grid*n_person), ones(1,n_grid))';
indx.iWeight2 = reshape(1:n_period*n_person, n_period,n_person);
indx.iWeight2 = kron(indx.iWeight2, ones(1,n_grid));
indx.iWeight2 = indx.iWeight2(:);
indx.jWeight2 = kron((1:n_grid*n_person), ones(1,n_period))';

DAT.ms = zeros(n_obs,n_mc);
DAT.ms0 = zeros(n_period,n_mc);

MAT.X = zeros(n_obs,n_theta1,n_mc);
MAT.Z = zeros(n_obs,n_inst,n_mc);
MAT.invPHI = zeros(n_inst,n_inst,n_mc);
MAT.XZPHIZXZ = zeros(n_theta1,n_obs,n_mc);
MAT.Xrc = zeros(n_obs,n_theta2,n_mc);

RND.nu = zeros(n_theta2,n_person,n_mc);
RND.seed = zeros(n_theta2,n_seed,n_mc);

end

