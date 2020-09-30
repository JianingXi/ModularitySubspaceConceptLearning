addpath('./bin_ROC');

network_name = 'STRING'; % or 'iRefIndex'
mutation_data_list = 'prostate_cancer_data'; % or 'thyroid_cancer_data'

% -- input: Data --
load(['./Network/' network_name '.mat']);
Data_dir = ['./Data/' mutation_data_list '_aligned_to_' network_name];
load([Data_dir '/mutation_data.mat']);

% -- Modularity Subspace Learning --
lambda_X = 1;
lambda_H = 100;
K_dim = 64;
lambda_eps = 100*lambda_H;

[W_mat,H_mat,H_eps,W_init] = ...
    f01_module_subspace(D_mut,Adj_mat,K_dim,lambda_X,lambda_H,lambda_eps,W_init);

% -- Bayesian concept learning --
load([Data_dir '/CV_indices.mat']);
N_cv = 10;

TPR = cell(N_cv,1);
FPR = cell(N_cv,1);
AUC_ROC = zeros(N_cv,1);
% -- Cross Validation --
for i_cv = 1:N_cv
    ind_test = [CV_cell{i_cv,1}; CV_cell{i_cv,2}];
    ind_train = setdiff(1:length(Adj_mat),ind_test)';

    % -- Training --
    Mdl = fitcnb(H_mat(ind_train,:),TrueIDX(ind_train,:),...
                'DistributionNames','kernel','Kernel','box');

    % -- Testing --
    [~, score_t] = predict(Mdl,H_mat(ind_test,:));
    Score_CV = score_t(:,2);

    % -- ROC evaluation --
    [TPR{i_cv},FPR{i_cv},AUC_ROC(i_cv)] = ...
                ROC_curve(TrueIDX(ind_test),Score_CV);
end
TPR_mean = CombineROC_CV(TPR);
FPR_mean = CombineROC_CV(FPR);
AUC_ROC_mean = mean(AUC_ROC);

plot(FPR_mean,TPR_mean,'Color','r');
grid on; box on; xlim([0 1]); ylim([0 1]);

rmpath('./bin_ROC');
        