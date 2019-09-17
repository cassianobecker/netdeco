function out = solve_discrete_tnn_feasibility(A0, B, par)
%#ok<*EQEFF>
%#ok<*VUNUS>

n    = size(A0,1);
I_n  = eye(n);
I_2n = eye(2*n);
BBT  = B*B';

% ----- initialization

max_abs_lam_A0 = max(abs(eig(A0)));

if max_abs_lam_A0 >1
     A0t = A0/(1e-3 + max_abs_lam_A0);
     W0 = dlyap(A0t, BBT);
else
     W0 = dlyap(A0, BBT);
end

if max_abs_lam_A0 >1
    fprintf('!!! Warning: Initial unstable A0 matrix, max(abs(eig(A0))))=%1.4e', max_abs_lam_A0)
end

%W0 = cov(rand(20,n))/100;
%W0 = zeros(n);

H0 = W0*A0';

Q  = [-BBT; zeros(n)];

M0 = [
     H0'  , -W0 ; ...
    -W0  ,   H0 ];

N0 = [ A0' ; I_n];

Z0 = [
       Q   ,    M0 ; ...
      N0   ,  I_2n ; ...
    ];

[U,~,V] = svd(Z0);

K1 = U(:,1:2*n)';
K2 = V(:,1:2*n)';

fprintf('\n##################### BEGIN OPTIMIZATION ###################\n');

% ----- start history

hst0.A0   = A0;
hst0.B    = B;

lams0 = eig(W0);

hst0.max_abs_lam_A = max(abs(eig(A0)));
hst0.is_schur = (hst0.max_abs_lam_A < 1.0);

hst0.ctrb = ctrb(A0, B);
hst0.rank_ctrb = rank(hst0.ctrb);
hst0.is_ctrb = (hst0.rank_ctrb == n);

hst0.min_lam   = min(lams0);
hst0.sum_lam_k = sum(lams0(1:par.s.k_lams));
hst0.tr_inv    = sum(1./lams0);

hst0.min_lam_bar   = par.s.min_lam_bar;
hst0.sum_lam_k_bar = par.s.sum_lam_k_bar;
hst0.tr_inv_bar    = par.s.tr_inv_bar;

gk = 1;

no_dec_ctr = 0;

fprintf('| iter %3d of %3d ', 0, par.m.MAX_ITER);

fprintf('| tnn: %+.3e', tnn(Z0, 2*n));

if par.s.mode.do_min_lam

    fprintf('| orig_lam_1: %+.3e| cur_lam_1: %+.3e| tar_lam_1: %+.3e', ...
        hst0.min_lam, hst0.min_lam, hst0.min_lam_bar);
end

if par.s.mode.do_tr_inv

    fprintf('| orig_tr_inv: %+.3e | cur_tr_inv: %+.3e | tar_tr_inv: %+.3e ', ...
        hst0.tr_inv , hst0.tr_inv , hst0.tr_inv_bar);
end

fprintf('|\n');


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
    
    hst(gk).max_abs_lam_A = max(abs(eig(A0 + opt.D)));
    hst(gk).is_schur = (hst(gk).max_abs_lam_A < 1.0);
    
    hst(gk).ctrb = ctrb(A0 + opt.D, B);
    hst(gk).rank_ctrb = rank(hst(gk).ctrb);
    hst(gk).is_ctrb = (hst(gk).rank_ctrb == n);
    
    if min(lams) < 0       
       fprintf('Rank Ctrb = %d of %d\n', hst(gk).rank_ctrb, size(A0,1));
    end

    hst(gk).lambda_1{k}  = min(lams);
    hst(gk).lambdas{k}   = lams;
    hst(gk).tr_inv{k}    = sum(1./lams);
    hst(gk).sum_lam_k{k} = sum(lams(1:par.s.k_lams));
    
    hst(gk).tnn{k}     = tnn(opt.Z, 2*n);
    hst(gk).status     = opt.status;
    hst(gk).optval     = opt.optval;
    
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
    
    fprintf('| status: %8s| opt_val: %+1.3e', opt.status, opt.optval);
    
    fprintf('|\n');
    
    % -------- Stopping Criteria ------------

    if par.s.mode.do_min_lam
        
        if hst(gk).lambda_1{k} >= hst0.min_lam_bar && hst(gk).tnn{k} < par.m.tol_tnn
            reason = 'SUCCESS: inequality residual below desired tolerance';
            fprintf('>>>> %s L\n', reason)
            break
        end

    end
    
    if par.s.mode.do_tr_inv
        
        if  hst(gk).tr_inv{k} <= hst0.tr_inv_bar && hst(gk).tnn{k} < par.m.tol_tnn
            reason = 'SUCCESS: inequality residual below desired tolerance';
            fprintf('>>>> %s L\n', reason)
            break
        end

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
            reason = 'LOW DESCREASE between iterations - STOPPED';
            fprintf('>>>> %s \n', reason)
            break
        end
    end
    
end

if k==par.m.MAX_ITER
   reason = sprintf('MAXIMUM number of iterations %d reached - STOPPED',k);
   fprintf('>>>> %s \n', reason)
end

%sum(sum(abs(D - par.s.amax) < par.m.tol_sparsity))
%sum(sum(abs(D - par.s.amin) < par.m.tol_sparsity))

%fprintf('Min. eval before/after: %1.3e / %1.3e\n', min(eig(W0)), min(eig(Wk)));

hst(gk).norm_D = norm(vec(opt.D),1); %#ok<*AGROW>
hst(gk).card_D = sum(abs(vec(opt.D)) > par.m.tol_sparsity);
hst(gk).D      = opt.D;
hst(gk).reason = reason;

fprintf('Norm of D: %1.3f, ', hst(gk).norm_D);
fprintf('Card of D: %1.3f\n', hst(gk).card_D);

fprintf('Min abs lam of A + D: %1.3f\n', hst(gk).max_abs_lam_A);
fprintf('Margin of stability of A + D: %1.3f\n', 1 - hst(gk).max_abs_lam_A);
fprintf('Controllability Rank of A + D, B: %1.3f\n', hst(gk).rank_ctrb);

out.hst  = hst;
out.hst0 = hst0;

end