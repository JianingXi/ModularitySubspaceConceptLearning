function [W_mat,H_mat,H_eps,W_init] = ...
    f01_module_subspace(D_mut,Adj_mat,K_dim,lambda_X,lambda_H,lambda_eps,W_init)

    % X: genes x samples
    [N_sample, P_node] = size(D_mut);
    % -- load network --
    [N_node, N_node2] = size(Adj_mat);
    if N_node ~= N_node2
        error('???');
    end

    degree_vec_raw = sum(Adj_mat,2);
    A_net = diag(degree_vec_raw.^(-0.5))*Adj_mat*diag(degree_vec_raw.^(-0.5));
    % -- load network done --

    eps_d = 10^-5;
    if ~exist('W_init','var')
        W_init = max(ones(N_sample,K_dim),eps_d); % init W
    end
    H_init = max(pinv(W_init)*D_mut,eps_d)'; % init Z
    H_eps_init = zeros(size(H_init));
    W_prev = W_init; H_prev = H_init; H_eps_prev = H_eps_init;

    res_H = 2; res_W = 2; res_H_eps = 2;
    res_th = 10^-3;  cnt = 0;
    while res_W >= res_th || res_H >= res_th || res_H_eps >= res_th
        H_hat = max(H_prev + H_eps_prev, eps_d);

        W_mat = W_prev.* ((D_mut*H_hat) + eps_d) ./ ...
            (W_prev*(H_hat'*H_hat) + eps_d);

        H_mat = H_prev.*(N_sample*P_node*K_dim*A_net*H_prev + lambda_X*(P_node^2)*K_dim*D_mut'*W_mat + eps_d) ./ ...
            (ones(P_node,1)*sum(H_prev,1)*(N_sample*K_dim/2) + lambda_X*(P_node^2)*K_dim*H_hat*(W_mat'*W_mat)...
            + lambda_H*N_sample*P_node^2*H_prev + eps_d);
        
        H_eps = ( (lambda_X*K_dim*(W_mat'*W_mat) + lambda_eps*N_sample*eye(K_dim)) \ ...
            (lambda_X*K_dim*(W_mat'*D_mut + (W_mat'*W_mat)*H_mat')) )';
        
        diag_H = sum(H_mat,1);
        diag_H(diag_H <= 0) = 1;
        H_mat = H_mat*diag(diag_H.^(-1));
        W_mat = W_mat*diag(diag_H);
        H_eps = H_eps*diag(diag_H);

        res_W = norm(W_mat-W_prev,'fro')/(norm(W_prev,'fro') + eps_d);
        res_H = norm(H_mat-H_prev,'fro')/(norm(H_prev,'fro') + eps_d);
        res_H_eps = norm(H_eps-H_eps_prev,'fro')/(norm(H_eps_prev,'fro') + eps_d);

        W_prev = W_mat; H_prev = H_mat; H_eps_prev = H_eps;
        
        cnt = cnt+1;
        if mod(cnt,100) == 0
            disp(cnt);
        end
        
    end % while

end % function
















