function run_ieeebus()

close all;

data_file    = 'ieeebus';
results_file = 'results_ieeebus';

fpath = fullfile(getRoot(),'data',data_file);
load(fpath)

if rank(ctrb(A0,B)) < size(A0, 1)
    error('system uncontrollable');
end

%%%%%%%%%%%%%%% ENUMERATE VARYING PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

par.s.min_lam_gains = [10, 100];
par.s.tr_inv_gains  = [1/10 1/50];

par.s.modes_tr_inv  = [0 1];
par.s.modes_min_lam = [1 0];

%%%%%%%%%%%%%%%% FIXED PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%

% structural contraints
par.s.amax =  0.5;
par.s.amin = -par.s.amax;

% penalization path D
par.s.base_gam = 10.0;
par.s.min_gam  = -3;
par.s.max_gam  = -1;
par.s.num_gam  = 40;

%%%%%%%%%%%%%%% SOLVER PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%

par.m.opt_mode = 'penalize_D';

par.m.tol_tnn  = 1e-7;
par.m.tol_eq   = 1e-8;
par.m.tol_sparsity = 1e-4;
par.m.rel_tol_dec = 1e-4;

par.m.max_no_decrease = 8;
par.m.MAX_ITER = 500;

%%%%%%%%%%%%%%%%% LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%

ni = length(par.s.min_lam_gains);
nj = length(par.s.modes_min_lam);

for i=1:ni
    
    for j=1:nj
        
        %%%%%%%%%%%%%%% SELECT PROBLEM DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        fprintf('\n\n');
        fprintf('#######################################################################################\n');
        fprintf('###################### SOLVING PROBLEM (%d, %d) of (%d, %d) ###############################\n', i, j, ni, nj);
        fprintf('#######################################################################################\n');
        
        %%%%%%%%%%%%%%% DEFINE PROBLEM PARAMTERS %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % target eigenvalues - MINIMUM LAMBDA
        par.s.mode.do_min_lam = par.s.modes_min_lam(i);
        par.s.gain_min_lam = par.s.min_lam_gains(j);
        
        % target eigenvalues - TRACE INVERSE
        par.s.mode.do_tr_inv = par.s.modes_tr_inv(i);
        par.s.gain_tr_inv = par.s.tr_inv_gains(j);
        
        % target eigenvalues - WORST K
        par.s.mode.do_sum_lam_k  = 0;
        par.s.k_lams = 3;
        par.s.gain_sum_lam_k = 5;
        
        %%%%%%%%%%%%%%%%% RUN OPTIMIZATION %%%%%%%%%%%%%%%%%
        
        out{i,j} = run_p_IEEEbus_unit(A0, B, par); %#ok<*AGROW>
        
        fig_reg_path(out{i,j});
        
    end
end

%%%%%%%%%%%%%%%% SAVE STATS %%%%%%%%%%%%%%%%%%%%%%%%%%%

fpath = fullfile(getRoot(), 'out', 'mat', results_file);

save(fpath, 'out');

end


function out = run_p_IEEEbus_unit(A0, B, par)

%%%%%%%%%%%%%%% TARGET PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%

W0 = dlyap(A0, B*B');
lams = sort(real(eig(W0)),'ascend');

par.s.min_lam = lams(1);
par.s.min_lam_bar = par.s.gain_min_lam*par.s.min_lam;

par.s.tr_inv = sum(1./lams);
par.s.tr_inv_bar   = par.s.gain_tr_inv*par.s.tr_inv;

par.s.sum_lam_k  = sum(lams(1:par.s.k_lams));
par.s.sum_lam_k_bar = par.s.gain_sum_lam_k*lams(par.s.k_lams);

%%%%%%%%%%%%%%% CALL SOLVER %%%%%%%%%%%%%%%%%%%%%%%%%%%

out = solve_discrete_tnn_penalty(A0, B, par);

out.par = par;

end