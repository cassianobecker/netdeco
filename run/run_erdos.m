function run_erdos()

close all;

problem_file = 'problems_erdos';
results_file = 'results_erdos_new';

fpath = fullfile(getRoot(), 'data', problem_file);
load(fpath)

%%%%%%%%%%%%%%% ENUMERATE VARYING PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

par.s.lam_min_gains = [10 50];
par.s.tr_inv_gains  = 1./[10 50];

par.s.modes_min_lam = [1 0];
par.s.modes_tr_inv  = [0 1];

nl = 100;
ni = length(par.s.modes_min_lam);
nj = length(par.s.lam_min_gains);

% aditional
par.s.k_lams = 3;
par.s.gain_sum_lam_k = 5;

% structural contraints
par.s.amax =  0.5;
par.s.amin = -par.s.amax;

%%%%%%%%%%%%%%% SOLVER PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%

par.m.opt_mode = 'feasibility';

par.m.tol_tnn  = 1e-7;
par.m.tol_eq   = 1e-7;
par.m.tol_sparsity = 1e-4;
par.m.rel_tol_dec = 1e-4;

par.m.max_no_decrease = 8;
par.m.MAX_ITER = 400;

for i=1:ni
    
    par.s.mode.do_min_lam    = par.s.modes_min_lam(i);
    par.s.mode.do_tr_inv     = par.s.modes_tr_inv(i);
    par.s.mode.do_sum_lam_k  = 0;
    
    for j=1:nj

        par.s.gain_min_lam = par.s.lam_min_gains(j);
        par.s.gain_tr_inv = par.s.tr_inv_gains(j);
        
        for l=1:nl
            
            %%%%%%%%%%%%%%% SELECT PROBLEM DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            fprintf('\n\n');
            fprintf('#######################################################################################\n');
            fprintf('#################### SOLVING PROBLEM (%d, %d, %04d) of (%d, %d, %04d) #####################\n', i , j, l, ni, nj, nl);
            fprintf('#######################################################################################\n');
            
            %%%%%%%%%%%%%%% DEFINE PROBLEM PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            offset = 0;
            A0 = As{offset+l};
            B  = Bs{offset+l};
            
            %%%%%%%%%%%%%%%%% RUN OPTIMIZATION %%%%%%%%%%%%%%%%%
            
            out{i, j, l} = run_p_erdos_unit(A0, B, par); %#ok<*AGROW>
            
        end
        
    end
    
    %%%%%%%%%%%%%%%% SAVE STATS %%%%%%%%%%%%%%%%%%%%%%%%%%%
end

fpath = fullfile(getRoot(), 'out', 'mat', results_file);
save(fpath, 'out');

end


function out = run_p_erdos_unit(A0, B, par)

%%%%%%%%%%%%%%% TARGET PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%

% unperturbed gramian spectrum
W0 = dlyap(A0, B*B');
lams = sort(real(eig(W0)),'ascend');

par.s.min_lam = lams(1);
par.s.min_lam_bar = par.s.gain_min_lam*par.s.min_lam;

par.s.tr_inv = sum(1./lams);
par.s.tr_inv_bar = par.s.gain_tr_inv*par.s.tr_inv;

par.s.sum_lam_k = sum(lams(1:par.s.k_lams));
par.s.sum_lam_k_bar = par.s.gain_sum_lam_k*lams(par.s.k_lams);

%%%%%%%%%%%%%%% CALL SOLVER %%%%%%%%%%%%%%%%%%%%%%%%%%%

out = solve_discrete_tnn_feasibility(A0, B, par);

out.par = par;

end