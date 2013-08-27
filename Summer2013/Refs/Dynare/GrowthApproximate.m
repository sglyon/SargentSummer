%
% Status : main Dynare file 
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

clear all
clear global
tic;
global M_ oo_ options_ ys0_ ex0_ estimation_info
options_ = [];
M_.fname = 'GrowthApproximate';
%
% Some global variables initialization
%
global_initialization;
diary off;
logname_ = 'GrowthApproximate.log';
if exist(logname_, 'file')
    delete(logname_)
end
diary(logname_)
M_.exo_names = 'e';
M_.exo_names_tex = 'e';
M_.endo_names = 'c';
M_.endo_names_tex = 'c';
M_.endo_names = char(M_.endo_names, 'k');
M_.endo_names_tex = char(M_.endo_names_tex, 'k');
M_.endo_names = char(M_.endo_names, 'lab');
M_.endo_names_tex = char(M_.endo_names_tex, 'lab');
M_.endo_names = char(M_.endo_names, 'z');
M_.endo_names_tex = char(M_.endo_names_tex, 'z');
M_.param_names = 'bet';
M_.param_names_tex = 'bet';
M_.param_names = char(M_.param_names, 'the');
M_.param_names_tex = char(M_.param_names_tex, 'the');
M_.param_names = char(M_.param_names, 'del');
M_.param_names_tex = char(M_.param_names_tex, 'del');
M_.param_names = char(M_.param_names, 'alp');
M_.param_names_tex = char(M_.param_names_tex, 'alp');
M_.param_names = char(M_.param_names, 'tau');
M_.param_names_tex = char(M_.param_names_tex, 'tau');
M_.param_names = char(M_.param_names, 'rho');
M_.param_names_tex = char(M_.param_names_tex, 'rho');
M_.param_names = char(M_.param_names, 's');
M_.param_names_tex = char(M_.param_names_tex, 's');
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 4;
M_.param_nbr = 7;
M_.orig_endo_nbr = 4;
M_.aux_vars = [];
M_.Sigma_e = zeros(1, 1);
M_.H = 0;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
erase_compiled_function('GrowthApproximate_static');
erase_compiled_function('GrowthApproximate_dynamic');
M_.lead_lag_incidence = [
 0 3 7;
 1 4 0;
 0 5 8;
 2 6 9;]';
M_.nstatic = 0;
M_.nfwrd   = 2;
M_.npred   = 1;
M_.nboth   = 1;
M_.equations_tags = {
};
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(4, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(7, 1);
M_.NNZDerivatives = zeros(3, 1);
M_.NNZDerivatives(1) = 18;
M_.NNZDerivatives(2) = 38;
M_.NNZDerivatives(3) = -1;
options_.periods = 1000;
M_.params( 1 ) = 0.987;
bet = M_.params( 1 );
M_.params( 2 ) = 0.357;
the = M_.params( 2 );
M_.params( 3 ) = 0.012;
del = M_.params( 3 );
M_.params( 4 ) = 0.4;
alp = M_.params( 4 );
M_.params( 5 ) = 2;
tau = M_.params( 5 );
M_.params( 6 ) = 0.95;
rho = M_.params( 6 );
M_.params( 7 ) = 0.007;
s = M_.params( 7 );
%
% INITVAL instructions
%
options_.initval_file = 0;
oo_.steady_state( 2 ) = 1;
oo_.steady_state( 1 ) = 1;
oo_.steady_state( 3 ) = 0.3;
oo_.steady_state( 4 ) = 0;
oo_.exo_steady_state( 1 ) = 0;
if M_.exo_nbr > 0;
	oo_.exo_simul = [ones(M_.maximum_lag,1)*oo_.exo_steady_state'];
end;
if M_.exo_det_nbr > 0;
	oo_.exo_det_simul = [ones(M_.maximum_lag,1)*oo_.exo_det_steady_state'];
end;
%
% SHOCKS instructions
%
make_ex_;
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (1)^2;
M_.sigma_e_is_diagonal = 1;
steady;
options_.irf = 40;
options_.periods = 1000;
var_list_=[];
info = stoch_simul(var_list_);
datasaver('simudata',[]);
save('GrowthApproximate_results.mat', 'oo_', 'M_', 'options_');


disp(['Total computing time : ' dynsec2hms(toc) ]);
disp([char(10) 'Dynare Preprocessor Warning(s) Encountered:']);
disp('     WARNING: GrowthApproximate.mod:3.1-13: periods: this command is now deprecated and may be removed in a future version of Dynare. Please use the ''periods'' option of the ''simul'' command instead.');
disp('     WARNING: GrowthApproximate.mod:40.13-21: dr_algo option is now deprecated, and may be removed in a future version of Dynare');
diary off
