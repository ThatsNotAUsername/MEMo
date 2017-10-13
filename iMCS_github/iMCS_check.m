function [ check ] = iMCS_check(iMCS, obj, network ,varargin )
% this function check that after knocking out reactions of iMCS the
% objective reaction cannot carry any flux. 

%INPUT : iMCS: each column represent an iMCS. A value of 1 means that
%reaction is present and a value of 0 means reaction is not in the IMCS.
%         obj: index of the objective reaction
%         network : metabolic network with the fields .S .rev .lb and .ub
%         optional:
%
%                   1) tolerance: for the given tolerance we decide if a reaction is active or
%                   not. Default is 1e-06
%                   2) solver: default is cplex. otherwise 'gurobi'


% we first build the MILP 

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


[n_mets n_rxns]=size(network.S);

[n_rxns n_iMCS]=size(iMCS);


%row vector of binary values 1 means the canditate is a cut set. else
%it's 0
check=zeros(1,n_iMCS);

rev=network.rev;

obj_vector=zeros(n_rxns,1);
obj_vector(obj)=ones(1,length(obj));

Model.A=sparse(double(network.S));
Model.obj=obj_vector
Model.lb = network.lb;
Model.ub = network.ub;
vtype=ones(1,n_rxns);
vtype(:)='C';


%%%%%%%%%%%%% GUROBI %%%%%%%%%%%%%%%%%%

if strcmp(solver,'gurobi') %use gurobi
    
    Model.sense= char(char('=')*ones(n_mets,1));
    Model.rhs=zeros(n_mets,1);
    Model.vtype=char(vtype);
    Model.modelsense='max';

    params.OutputFlag = 0; % gurobi is silent
    params.FeasibilityTol = tol;
    %params.IntFeasTol = tol;
    params.OptimalityTol = tol;
    params.MarkowitzTol = 0.5;

 
%for each reaction in a IMCS we have to put a 0 flux. If
%the reaction is not part of the IMCS but is irreversible we have the
%inequality vi>0

initial_A=Model.A;
initial_sense=Model.sense;
initial_rhs=Model.rhs;

    for i=1:n_iMCS
        
        cut_set=iMCS(:,i);
    
        for j=1:n_rxns
            if rev(j) %the reaction j is reversible 
            else vect=zeros(1,n_rxns); %the reaction is irreversible
                vect(j)=1; 
                Model.A(end+1,:)=vect;
                Model.rhs(end+1,:)= 0;
                if  cut_set(j) % the reaction is irreversible and is part of the iMCS
                     Model.sense(end+1,:)=char('=');
                else Model.sense(end+1,:)=char('>'); % the reaction is irreversible but is not part of the iMCS
                end
            end
        end
    
    %solve the MILP
    result = gurobi(Model,params);
    
     if strcmp(result.status, 'INFEASIBLE') 
          disp(' ');
            disp('#######################################');
            disp(' ');
            disp('Model is infeasible');
            disp(' ');
            disp('#######################################');
            disp(' ');
            check(i)=1;
     else res=result.objval
         if abs(res)<tol % the solution is considered to be equal to 0
             check(i)=1; 
         end
     end
     
Model.A=initial_A;
Model.sense=initial_sense;
Model.rhs=initial_rhs;

    end
   
 
%%%%%%%%%%%%%%%%%%%% CPLEX %%%%%%%%%%%%%%%%%

elseif strcmp(solver,'cplex')
    tosolve = Cplex;
    tosolve.Model = Model;
    tosolve.Model.lhs = zeros(n_mets,1);
    tosolve.Model.rhs = zeros(n_mets,1);
    tosolve.Model.ctype = char(vtype);
    tosolve.Model.sense='maximize'
    
    tosolve.Param.simplex.tolerances.optimality.Cur = tol; % default: 1e-06
    tosolve.Param.simplex.tolerances.feasibility.Cur = tol; % default: 1e-06
    %tosolve.Param.mip.tolerances.integrality.Cur = tol; % default: 1e-06
    tosolve.Param.mip.display.Cur = 'off'; % Solver should not show its iterations
    %tosolve.writeModel('cplex_Each_Ind.lp');
    tosolve.Param.emphasis.numerical.Cur = true;
    
    %for each reaction in a iMCS we have to put a 0 flux. If
    %the reaction is not part of the iMCS but is irreversible we have the
    %inequality vi>0
    

    
    initial_A=tosolve.Model.A;
    initial_rhs=tosolve.Model.rhs;
    initial_lhs=tosolve.Model.lhs;
    
    for i=1:n_iMCS
        cut_set=iMCS(:,i);
    
        for j=1:n_rxns
            if rev(j) %reaction is reversible 
            else vect=zeros(1,n_rxns); %reaction is irreversible
                vect(j)=1; 
                tosolve.Model.A(end+1,:)=vect;
                tosolve.Model.lhs(end+1,:)= 0;
                if  cut_set(j) % the reaction is irreversible and is part of an iMCS
                     tosolve.Model.rhs(end+1,:)=0;
                else tosolve.Model.rhs(end+1,:)=inf; % the reaction is irreversible but is not part of an iMCS
                end
            end
        end
        
        tosolve.solve();
        result = tosolve.Solution;
    
        if strcmp(result.status,'INFEASIBLE') | result.status ==103
            disp(' ');
            disp('#######################################');
            disp(' ');
            disp('Model is infeasible');
            disp(' ');
            disp('#######################################');
            disp(' ');
            check(i)=1;
        
        else res=result.objval;
            if abs(res)<tol %the solution is considered to be equal to 0
               check(i)=1;
            else fprintf('For iMCS %d the objective value is not 0 : %f \n',i,res)
            end
        end
        
    tosolve.Model.A=initial_A;
    tosolve.Model.rhs=initial_rhs;
    tosolve.Model.lhs=initial_lhs;
       
    end
    
   
else error('Wrong solver. Either cplex or gurobi.'); 
end    


if sum(check)==n_iMCS
    disp('')
    disp('All the sets are really IMCS')
    disp('')
else fprintf('There are %d sets that are not IMCS \n',n_iMCS-sum(check))
end
    
end

