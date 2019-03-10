function opt = cvx_tnn_feasibility(A0, B, K1, K2, par)
%#ok<*EQEFF>
%#ok<*STOUT>
%#ok<*VUNUS>

n    = size(A0,1);
I_n  = eye(n);
I_2n = eye(2*n);
S    = (abs(A0)==0);
Q    = [-B*B'; zeros(n)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning off

cvx_precision default

cvx_begin sdp quiet

variable M(2*n,2*n)
variable N(2*n,n)
variable Z(4*n,3*n)
variable W(n,n) semidefinite;
variable A(n,n);
variable D(n,n);
variable H(n,n);

if par.s.mode.do_tr_inv
    variable L2(n,n) semidefinite;
end

if par.s.mode.do_sum_lam_k
    variable s
    variable L1(n,n) semidefinite;
end

minimize(norm_nuc(Z) - trace(K1*Z*K2'))

Z == [
    Q , M    ; ...
    N , I_2n ; ...
    ];

M == [
    H'  , -W ; ...
    -W  ,  H
    ];

N == [ A' ; I_n ];

A == A0 + D;

if par.s.mode.do_min_lam
    
    W - par.s.min_lam_bar*I_n >= 0;
    
end

if par.s.mode.do_sum_lam_k
    
    - par.s.sum_lam_k_bar - par.s.k_lams*s - trace(L1) >= 0;
    
    L1 >= 0;
    
    L1 + W + s*I_n >= 0;
    
end

if par.s.mode.do_tr_inv
    
    par.s.tr_inv_bar - trace(L2) >= 0;
    
    [W I_n; I_n L2] >= 0;
    
end

vec(D)'*diag(vec(S)) == zeros(1, n^2);
vec(D) <=  par.s.amax*ones(n^2, 1);
vec(D) >=  par.s.amin*ones(n^2, 1);

cvx_end

opt.Z = Z;
opt.W = W;
opt.D = D;
opt.M = M;
opt.N = N;
opt.A = A;
opt.status = lower(cvx_status);
opt.optval = cvx_optval;

%%%%%%%%

end