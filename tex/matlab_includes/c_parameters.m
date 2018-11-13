%% Divide the data sets into overlapping bits 
set{1} = 1:10; 
set{2} = 10:20;
set{3} = 20:30; 
set{4} = 30:40; 
set{5} = 40:50;

% Choose the reconstruction parameters 
w{1} = 0.2 * ones(numel(set{1}),1);
w{2} = 0.3 * ones(numel(set{2}),1);
w{3} = 0.4 * ones(numel(set{3}),1);
w{4} = 0.3 * ones(numel(set{4}),1);
w{5} = 0.2 * ones(numel(set{5}),1);

alpha{1} = 500 * ones(numel(set{1}),1);
alpha{2} = 500 * ones(numel(set{2}),1);
alpha{3} = 500 * ones(numel(set{3}),1);
alpha{4} = 500 * ones(numel(set{4}),1);
alpha{5} = 500 * ones(numel(set{5}),1);

gamma{1} = 75 * ones(numel(set{1}),1);
gamma{2} = 75 * ones(numel(set{2}),1);
gamma{3} = 75 * ones(numel(set{3}),1);
gamma{4} = 75 * ones(numel(set{4}),1);
gamma{5} = 75 * ones(numel(set{5}),1);