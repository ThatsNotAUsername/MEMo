function [ hits,hits_card ] = hitting_sets( generators,reaction,cardinality,varargin)

% given a set of generators and one or several objective reaction we will find the
% corresponding hitting sets
% if reaction=0 we will search for hitting sets for the whole network

% we will solve this problem using a MILP:
% a is a binary vector with a_i = 1 if reaction is in hitting set, 0 else
% min sum a_i
% s.t. g^j * a >=1, where g^j is the j-th generator, which is a  given binary
% vector. >=1 makes sure that for each generator a reaciton is chosen
% no good cuts s.t. we compute in each step a different hitting set 


%INPUT : * generators:matrix of generator where one row represent one
%generator and one column represents one reaction. All the values are binaries and a reaction is active in a
%generator if it's set to 1.
%         * reaction : index of the reactions for which we search the
%         hitting sets
%         * cardinality : max cardinality of hitting set. If 0 we find all
%         the hitting sets.
%         * optional:
%
%                   1) tolerance: for the given tolerance we decide if a reaction is active or
%                   not. Default is 1e-06
%                   2) solver: default is cplex. otherwise 'gurobi'
%OUTPUT : * hits : a matrix where one column represent one imcs. A value of 1 means that the reaction is active 
%         * hits_card :  a vector of size cardinality where hits_card(j) is the number of hits of cardinality j

if nargin >= 4
    tol = varargin{1};
    if nargin == 5
        solver = varargin{2};
    else
        solver = 'cplex';
    end
else
    tol = 1e-06;
    solver = 'cplex';
end


[n_generators n_rxns] = size(generators);
use_generators_vector = zeros(1,n_generators);
if reaction % one or several objective reactions are given
    for i=1:n_generators
        g = generators(i,:);
        if sum(g(reaction))>0 %this means that at least one of the objective reactions is contained in g
            use_generators_vector(i) = 1;
        end
    end
else
    use_generators_vector(:) = 1;
end

use_generators = generators(logical(use_generators_vector),:);%use_generators are the generators containig the objective reactions

n_use_generators = nnz(use_generators_vector);

%fprintf('the number of generators used is %d \n',n_use_generators)

if not(n_use_generators)
    error('There are no generrators used. Variable n_use_generators is equal to 0 in hitting_sets.m')
end

%build the MILP

Model.A=sparse(double(use_generators));
Model.obj=ones(n_rxns,1);
Model.lb = zeros(n_rxns,1);
Model.ub = ones(n_rxns,1);
vtype=ones(1,n_rxns);
vtype(:)='B';

hits=[];
hits_card=zeros(1,cardinality);

found_all=1;

n_sol=0;

%%%%%%%%%%%%%%% GUROBI %%%%%%%%%%%%%%%%%%%%

if strcmp(solver,'gurobi') %use gurobi
    Model.sense= char(char('>')*ones(n_use_generators,1));
    Model.rhs=ones(n_use_generators,1);
    Model.vtype=char(vtype);
    Model.modelsense = 'min';
    
    params.OutputFlag = 0; % gurobi is silent
    params.FeasibilityTol = tol;
    %params.IntFeasTol = tol;
    params.OptimalityTol = tol;
    params.MarkowitzTol = 0.5;

    %Found the solution. For each set found we add a no good cut.
    
    while found_all

        result = gurobi(Model,params);

        if strcmp(result.status, 'INFEASIBLE') 
            disp(' ');
            disp('#######################################');
            disp(' ');
            disp('Found all IMCSs');
            disp('#######################################');
            disp(' ');
            found_all=0;
            break;
    
        else  cut_set=result.x;
            
            computed_card=nnz(result.x);
        
            if cardinality < computed_card
                disp(' ');
                disp('#######################################');
                disp(' ');
                disp('Found all IMCSs of given cardinality');
                disp(' ');
                disp('#######################################');
                disp(' ');
                found_all=0;
                break;
        
       
            else hits=[hits cut_set];
                hits_card(computed_card)= hits_card(computed_card)+1;
                %add constraint
                Model.A=[Model.A; transpose(cut_set)];
                Model.rhs=[Model.rhs;computed_card-1];
                Model.sense(end+1,1)='<';
             
                n_sol=n_sol+1;
            
            end
        
        end
    end
    
%%%%%%%%%%%%%%%% CPLEX %%%%%%%%%%%%%%%%%%%%

elseif strcmp(solver,'cplex') %use cplex
    tosolve = Cplex;
    tosolve.Model = Model;
    tosolve.Model.rhs = inf*ones(n_use_generators,1);
    tosolve.Model.lhs = ones(n_use_generators,1);
    tosolve.Model.ctype = char(vtype);
    tosolve.Model.sense='minimize';
    
    tosolve.Param.simplex.tolerances.optimality.Cur = tol; % default: 1e-06
    tosolve.Param.simplex.tolerances.feasibility.Cur = tol; % default: 1e-06
    %tosolve.Param.mip.tolerances.integrality.Cur = tol; % default: 1e-06
    tosolve.Param.mip.display.Cur = 'off'; % Solver should not show its iterations
    %tosolve.writeModel('cplex_Each_Ind.lp');
    tosolve.Param.emphasis.numerical.Cur = true;
    
    
    
    while found_all
    tosolve.solve();
   
    result = tosolve.Solution;
    
    if strcmp(result.status,'INFEASIBLE') | result.status ==103
        disp(' ');
        disp('#######################################');
        disp(' ');
        disp('Found all iMCSs ');
        disp(' ');
        disp('#######################################');
        disp(' ');
        found_all = 0;
        break;
    end
    cut_set = result.x;
    computed_card = nnz(cut_set);
    if cardinality < computed_card
        disp(' ');
        disp('#######################################');
        disp(' ');
        disp('Found all iMCSs of given cardinality');
        disp(' ');
        disp('#######################################');
        disp(' ');
        found_all = 0;
        break;
    else
    hits = [hits  cut_set];
    hits_card(computed_card)=hits_card(computed_card)+1;
    tosolve.Model.A(end+1,:) = logical(cut_set)';
    tosolve.Model.rhs(end+1,:) = computed_card -1 ;
    tosolve.Model.lhs(end+1,:) = -inf;
    n_sol=n_sol+1;
    end
    end
    
    
    
else error('Wrong solver. Either cplex or gurobi.');
    
end
    

disp(' ');
fprintf('The number of solution is %d',n_sol);

disp(' ');

end

