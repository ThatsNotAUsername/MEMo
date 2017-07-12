function [blocked_rxns, rev] = correctModel(network,varargin)
% Correct the reversibilities and find blocked reactions of the network.

% for that we will loop over the reactions and do an lp to check if the
% reactions are reversible or blocked.
clear tic
clear toc
tic
lb = network.lb;
ub = network.ub;

S = network.S;
if nargin == 2
    tol = varargin{1};
else
    tol = 1e-06;
end
[n_mets n_rxns] = size(S);

rev = lb<0; % whenever the lowerbound is less than zero the rxn i is potentially reversible, e.g. in rev is at this position a 1
blocked_rxns = zeros(n_rxns,1);

% use cplex

tosolve = Cplex;
tosolve.Model.obj = zeros(1,n_rxns);
tosolve.Model.lb = lb;
tosolve.Model.ub = ub;
tosolve.Model.A = sparse(S);
tosolve.Model.lhs = zeros(n_mets,1);
tosolve.Model.rhs = zeros(n_mets,1);
vtype = ones(1,n_rxns);
vtype(:) = 'C';
tosolve.Model.ctype = char(vtype);

tosolve.Param.simplex.tolerances.optimality.Cur = tol; % default: 1e-06
%tosolve.Param.simplex.tolerances.feasibility.Cur = 1e-09; % default: 1e-06
tosolve.Param.mip.tolerances.integrality.Cur = tol; % default: 1e-06
tosolve.Param.mip.display.Cur = 'off'; % Solver should not show its iterations
%tosolve.writeModel('cplex_Each_Ind.lp');
tosolve.Param.emphasis.numerical.Cur = true;

for ind_rxn = 1:n_rxns
    if rev(ind_rxn) % if ind_rxn is potentially reversible we minimize the flux and check if it really can have negative fluxvalue
        tosolve.Model.obj = zeros(1,n_rxns);
        tosolve.Model.obj(ind_rxn) = 1;
        tosolve.Model.sense = 'minimize';
        tosolve.solve();
        result = tosolve.Solution;
        obj_val = result.objval;
       if obj_val>=0 % if it is not reversible
           rev(ind_rxn) = 0;
           tosolve.Model.sense = 'maximize'; % check if it is blocked
           tosolve.solve();
           result = tosolve.Solution;
           obj_val = result.objval;
           if obj_val == 0
               blocked_rxns(ind_rxn) = 1;
           end
       end
    else
        % if the reaction is not reversible we check if it is blocked:
        % for that we first maximize and if the max is zero we than
        % minimize
        tosolve.Model.obj = zeros(1,n_rxns);
        tosolve.Model.obj(ind_rxn) = 1;
         n_nonZeros = nnz(tosolve.Model.obj);
        tosolve.Model.sense = 'maximize';
        tosolve.solve();
        result = tosolve.Solution;
        obj_val = result.objval;
        if obj_val == 0
            tosolve.Model.sense = 'minimize';
            tosolve.solve();
            result = tosolve.Solution;
            obj_val = result.objval;
            if obj_val == 0
                blocked_rxns(ind_rxn) = 1;
            end
        end
    end
end
       
time = toc;
disp(['time for correctReversibilities: ', num2str(time)]);
    
