function param = paramfile_tv

param           = struct;
param.norm_type = 'iso';
param.pdtol     = 1e-5;
param.niter     = 5e4;
param.int       = 100;

param.tau       = 0.3; 
param.sigma     = 0.3;

param.show      = true;

end