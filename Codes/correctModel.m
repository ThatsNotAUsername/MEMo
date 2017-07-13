function [new_network, blocked_rxns, change, blocked_mets] = correctModel(network,varargin)
% Correct the reversibilities and find blocked reactions of the network.
% if a reaction can only carry flux in a negative direction we replace hte
% corresponding column in the stoichiometric matrix by its negative.
% we also look for dead-end metabolites, metabolites which are not used, and remove them 

% for that we will loop over the reactions and do an lp to check if the
% reactions are reversible or blocked. The metabolites which are not used,
% zero-ros, will be removed at the end

%% Input:
% network: struct, with:
% -S:               stoichiometric matrix
% -lb:              lowerbounds for the reactionrates
% -ub:              upperbounds for the reactionrates
%
% optional:
%
% 1) tolerance:     for the given tolerance we decide if a reaction is active or
%                   not. Default is 1e-06
% 2) solver:        default is cplex. otherwise 'gurobi' or 'glpk'

%% Output
%-new_model: Struct as before. 
%   -rev            is the updated reversible vector
%   - S:            is the stoichiometric matrix, where the columns for 
%                   reactions which can carry only negative flux in the 
%                   original model are replaced by its negative. The    
%                   blocked reactions were removed
%   - lb:           same as before, but the blocked reactions were removed. 
%                   If the i-th reaction is not reversible, lb(i) is set to 0.
%   - ub:           same as before, but the blocked reactions were removed. 
%   - rxns:         blocked reactions are removed
%   - rxnNames:     blocked reactions are removed
%   - mets:         blocked metabolites are removed
%   - metNames:     blocked metabolites are removed
%- blocked_rxns:    0-1 vector, with 1 indicating that the i-th reaction is
%                   blocked
% change:           0-1 vector. 1 is indicating that the direction of the 
%                   ith reaction was change in the stoichiometric matrix. 

clear tic
clear toc
tic
lb = network.lb;
ub = network.ub;

S = network.S;
if nargin >= 2
    tol = varargin{1};
    if nargin == 3
        solver = varargin{2};
    else
        solver = 'cplex';
    end
else
    tol = 1e-06;
end

[n_mets n_rxns] = size(S);

rev = lb<0; % whenever the lowerbound is less than zero the rxn i is potentially reversible, e.g. in rev(i) = 1
blocked_rxns = zeros(n_rxns,1); % 1, if reactions is blocked.
irr = zeros(n_rxns,1); % 1, if reaction is irreversible
change = zeros(n_rxns,1); % 1, if reaction is irreversible, but can only carry flux in the reverse direction. We will replace the corresponding column of the stoichiometric matrix by its negative


Model.obj = zeros(1,n_rxns);
Model.lb = lb;
Model.ub = ub;
Model.A = sparse(S);
Model.rhs = zeros(n_mets,1);
vtype = ones(1,n_rxns);
vtype(:) = 'C';

if strcmp(solver,'cplex')
    % use cplex
    tosolve = Cplex;
    tosolve.Model = Model;
    tosolve.Model.lhs = zeros(n_mets,1);
    %tosolve.Model.rhs = zeros(n_mets,1);
    tosolve.Model.ctype = char(vtype);
    
    tosolve.Param.simplex.tolerances.optimality.Cur = tol; % default: 1e-06
    tosolve.Param.simplex.tolerances.feasibility.Cur = tol; % default: 1e-06
    %tosolve.Param.mip.tolerances.integrality.Cur = tol; % default: 1e-06
    tosolve.Param.mip.display.Cur = 'off'; % Solver should not show its iterations
    %tosolve.writeModel('cplex_Each_Ind.lp');
    tosolve.Param.emphasis.numerical.Cur = true;
    
elseif strcmp(solver,'gurobi')
    Model.sense = '=';
    Model.vtype = char(vtype);
    params.outputflag = 0; % gurobi is silent
    params.FeasibilityTol = tol;
    %params.IntFeasTol = tol;
    params.OptimalityTol = tol;
    params.MarkowitzTol = 0.5;
elseif strcmp(solver,'glpk')
    ctype = zeros(n_mets,1);
    ctype(:)='S';
 %   Model.vartype = char(vtype);
    params.msglev = 0; % glpk is silent
    params.tolbnd = tol;
    params.tolobj = tol;
else
    error('Wrong solver. Either cplex, gurobi or glpk.');
end

for ind_rxn = 1:n_rxns
    if rev(ind_rxn) % if ind_rxn is potentially reversible we minimize the flux and check if it really can have negative fluxvalue
        if strcmp(solver,'cplex')
            tosolve.Model.obj = zeros(1,n_rxns);
            tosolve.Model.obj(ind_rxn) = 1;
            tosolve.Model.sense = 'minimize';
            tosolve.solve();
            result = tosolve.Solution;
        elseif strcmp(solver,'gurobi')
            Model.obj = zeros(1,n_rxns);
            Model.obj(ind_rxn) = 1;
            Model.modelsense = 'min';
            result = gurobi(Model,params);
        elseif strcmp(solver,'glpk')
            Model.obj = zeros(1,n_rxns);
            Model.obj(ind_rxn) = 1;
            [xopt, result.objval] = glpk(Model.obj, Model.A , Model.rhs, Model.lb, Model.ub, char(ctype), char(vtype), 1, params);
        end
        min_obj_val = result.objval;
        if strcmp(solver,'cplex')
            tosolve.Model.sense = 'maximize'; % check if it is blocked
            tosolve.solve();
            result = tosolve.Solution;
        elseif strcmp(solver,'gurobi')
            Model.modelsense = 'max';
            result = gurobi(Model,params);
        elseif strcmp(solver,'glpk')
            [xopt, result.objval] = glpk(Model.obj, Model.A , Model.rhs, Model.lb, Model.ub, char(ctype), char(vtype), -1, params);
        end
        max_obj_val = result.objval;
        if min_obj_val<= -tol
            if max_obj_val >= tol % reaction is reversible
                rev(ind_rxn) = 1;
            else                  % reaction is irreversible but the direction is a different one
                change(ind_rxn) = 1;
                irr(ind_rxn) = 1;
                ub(ind_rxn) = -lb(ind_rxn);
                lb(ind_rxn) = 0;
            end
        elseif min_obj_val >= -tol
            if max_obj_val >= tol % reaction is irreversible
                irr(ind_rxn) = 1;
                lb(ind_rxn) = 0;
            else                  % reaction is blocked
                lb(ind_rxn) = 0;
                ub(ind_rxn) = 0;
                blocked_rxns(ind_rxn) = 1;
            end
        end       
    else
        % if the reaction is not reversible we check if it is blocked:
        % for that we first maximize and if the max is zero we than
        % minimize
        if strcmp(solver,'cplex')
            tosolve.Model.obj = zeros(1,n_rxns);
            tosolve.Model.obj(ind_rxn) = 1;
            tosolve.Model.sense = 'maximize';
            tosolve.solve();
            result = tosolve.Solution;
        elseif strcmp(solver,'gurobi')
            Model.obj = zeros(1,n_rxns);
            Model.obj(ind_rxn) = 1;
            Model.modelsense = 'max';
            result = gurobi(Model,params);
        elseif strcmp(solver,'glpk')
            Model.obj = zeros(1,n_rxns);
            Model.obj(ind_rxn) = 1;
            [xopt, result.objval] = glpk(Model.obj, Model.A , Model.rhs, Model.lb, Model.ub, char(ctype), char(vtype), -1, params);
        end
        obj_val = result.objval;
        if obj_val >= tol
            irr(ind_rxn) = 1;
        else
            blocked_rxns(ind_rxn) = 1;
            ub(ind_rxn) = 0;
        end
    end
end

% use all the gathered information:
S(:, logical(change)) = -S(:, logical(change)); % replace the columns of the changed reactions by the negative column
S(:,logical(blocked_rxns)) = [];
lb(logical(blocked_rxns)) = [];
ub(logical(blocked_rxns)) = [];
rev(logical(blocked_rxns)) = [];
new_network = network;
new_network.S = S;
new_network.lb = lb;
new_network.ub = ub;
new_network.rev = rev;
%try % if rxns exist
    new_network.rxns(logical(blocked_rxns)) = [];
%catch
%end
try % if rxnNames exist
    new_network.rxnNames(logical(blocked_rxns)) = [];
catch
end
% check for zero rows in S:
blocked_mets = ~any(S,2); 
new_network.S(logical(blocked_mets),:) = [];
try % if mets exist
    new_network.mets(logical(blocked_mets)) = [];
catch
end
try % if metNames exist
    new_network.metNames(logical(blocked_mets)) = [];
catch
end
time = toc;
disp(['time for correctReversibilities: ', num2str(time)]);
    
