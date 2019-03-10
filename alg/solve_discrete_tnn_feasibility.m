function out = solve_discrete_tnn_feasibility(A0, B, par)
%#ok<*EQEFF>
%#ok<*VUNUS>

n    = size(A0,1);
I_n  = eye(n);
I_2n = eye(2*n);
BBT  = B*B';

% ----- initialization

W0     = dlyap(A0, BBT);

H0 = W0*A0';

Q  = [-BBT; zeros(n)];

M0 = [
    H0'  , -W0 ; ...
    -W0  ,  H0 ];

N0 = [ A0' ; I_n];

Z0 = [
    Q   , M0 ; ...
    N0  , I_2n ; ...
    ];

[U,~,V] = svd(Z0);

K1 = U(:,1:2*n)';
K2 = V(:,1:2*n)';

fprintf('\n##################### BEGIN OPTIMIZATION ###################\n');

% ----- start history

hst0.A0   = A0;
hst0.B    = B;

lams0 = eig(W0);

hst0.min_lam   = min(lams0);
hst0.sum_lam_k = sum(lams0(1:par.s.k_lams));
hst0.tr_inv    = sum(1./lams0);

hst0.min_lam_bar   = par.s.min_lam_bar;
hst0.sum_lam_k_bar = par.s.sum_lam_k_bar;
hst0.tr_inv_bar    = par.s.tr_inv_bar;

gk = 1;

no_dec_ctr = 0;

for k=1:par.m.MAX_ITER
    
    fprintf('| iter %3d of %3d ', k, par.m.MAX_ITER);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    opt = cvx_tnn_feasibility(A0, B, K1, K2, par);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [U,~,V] = svd(full(opt.Z));
    K1 = U(:, 1:2*n)';
    K2 = V(:, 1:2*n)';
    
    % ------- history ------------
    
    Wk = dlyap(A0 + opt.D, BBT);
    lams = sort(eig(Wk), 'ascend');
    
    hst(gk).W{k}        = Wk;
    hst(gk).W_var{k}    = opt.W;
    
    hst(gk).lambda_1{k}  = min(lams);
    
    if min(lams) < 0
       CAB = ctrb(A0 + opt.D, B);
       fprintf('Rank Ctrb = %d of %d\n', rank(CAB), size(A0,1));
    end
    
    hst(gk).lambdas{k}   = lams;
    hst(gk).tr_inv{k}    = sum(1./lams);
    hst(gk).sum_lam_k{k} = sum(lams(1:par.s.k_lams));
    
    hst(gk).tnn{k}     = tnn(opt.Z, 2*n);
    
    fprintf('| status: %8s| opt_val: %+1.3e', opt.status, opt.optval);
    fprintf('| tnn: %+.3e', hst(gk).tnn{k});
    
    if par.s.mode.do_min_lam
        
        fprintf('| orig_lam_1: %+.3e| cur_lam_1: %+.3e| tar_lam_1: %+.3e', ...
            hst0.min_lam, hst(gk).lambda_1{k}, hst0.min_lam_bar);
    end

    if par.s.mode.do_sum_lam_k
        
        fprintf('| orig_slam_k: %+.3e | cur_slam_k: %+.3e | tar_slam_k: %+.3e ', ...
            hst0.sum_lam_k, hst(gk).sum_lam_k{k}, hst0.sum_lam_k_bar);
    end
    
    if par.s.mode.do_tr_inv
        
        fprintf('| orig_tr_inv: %+.3e | cur_tr_inv: %+.3e | tar_tr_inv: %+.3e ', ...
            hst0.tr_inv, hst(gk).tr_inv{k}, hst0.tr_inv_bar);
    end
    
    fprintf('|\n');
    
    % -------- Stopping Criteria ------------

    if hst(gk).tnn{k} < par.m.tol_tnn
        fprintf('# SUCCESS: inequality residual below desired tolerance\n')
        break
    end
    
    if k > 1
        
        if (hst(gk).tnn{k-1} - hst(gk).tnn{k}) / hst(gk).tnn{k-1} <= par.m.rel_tol_dec
            no_dec_ctr = no_dec_ctr + 1;
            fprintf('Setting NO_DECREASE counter to %d of %d\n', no_dec_ctr, par.m.max_no_decrease);
        elseif no_dec_ctr > 0
            no_dec_ctr = 0;
            fprintf('Resetting NO_DECREASE counter\n');
        end
        
        if no_dec_ctr >= par.m.max_no_decrease
            fprintf('>>>> LOW DESCREASE between iterations - STOPPED\n')
            break
        end
    end
    
end

hst(gk).norm_D = norm(vec(opt.D),1); %#ok<*AGROW>
hst(gk).card_D = sum(abs(vec(opt.D)) > par.m.tol_sparsity);
hst(gk).D      = opt.D;

fprintf('Norm of D: %1.3f, ', hst(gk).norm_D);
fprintf('Card of D: %1.3f\n', hst(gk).card_D);

out.hst  = hst;
out.hst0 = hst0;

end