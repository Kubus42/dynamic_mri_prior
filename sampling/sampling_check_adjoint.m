%% Check, whether the operators sampling_op.m and sampling_adj.m are adjoint
clear; clc;

img = randn(128,256) + 1i * randn(128,256);
geom = sampling_geom(img,'spokes5','show',0);

err = 0; 
for i= 1:500
    data = sampling_op(img,geom);
    y = randn(size(data)) + 1i * randn(size(data)); 
    x = sampling_adj(y,geom);
    err = err + trace(img * x') - trace(data * y');
end

fprintf('Error: %6.6d.\n' , err);






