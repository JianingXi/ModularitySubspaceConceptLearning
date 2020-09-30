

function Mean_out = CombineROC_CV(res_score_cur)
    N_cv = length(res_score_cur);
    
    len_vec = cellfun(@length,res_score_cur);
    max_len = max(len_vec);
    mat_cur = ones(max_len,N_cv);
    
    for i_cv = 1:N_cv
        ind_in = 1:len_vec(i_cv);
        mat_cur(ind_in,i_cv) = res_score_cur{i_cv};
    end
    Mean_out = mean(mat_cur,2);
end





% 
% function [Mean_out, Std_out] = CombineROC_CV(res_score_cur)
%     N_cv = length(res_score_cur);
%     
%     if N_cv < 2
%         Mean_out = res_score_cur{1};
%     elseif N_cv == 2
%         Mean_out = CombineTwoTPRVec(res_score_cur{1},res_score_cur{2});
%     elseif N_cv > 2
%         Mean_out = CombineTwoTPRVec(res_score_cur{1},res_score_cur{2});
%         for i_cv = 3:N_cv
%             Mean_out = CombineTwoTPRVec(res_score_cur{i_cv},Mean_out);
%         end
%     end
%     
%     len_final = length(Mean_out);
%     mat_cur = zeros(len_final,N_cv);
%     for i_cv = 1:N_cv
%         mat_cur(:,i_cv) = 2*CombineTwoTPRVec(res_score_cur{i_cv},zeros(len_final,1));
%     end
%     Std_out = std(mat_cur,[],2);
% end
% 
% function y_out = CombineTwoTPRVec(f1,f2)
%     % [0, 1]
% 
%     L1 = length(f1);
%     L2 = length(f2);
%     LN = lcm(L1,L2);
% 
%     y1 = interp1q((0:(L1-1))'/(L1-1),f1,(0:(LN-1))'/(LN-1)); 
%     y2 = interp1q((0:(L2-1))'/(L2-1),f2,(0:(LN-1))'/(LN-1)); 
%     y_out = 0.5*(y1+y2);
% end