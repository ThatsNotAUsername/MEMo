function [ iMCS, iMCS_card ] = iMCS_computation( network,MMB,cardinality,tol,solver,varargin )
%Compute the irreversible minimal cut sets 
%INPUT : network : a structure with the fields .rxnNames .S .rev .lb .ub 
%        MMB : A matrix where one row is one mmb. A value of 1 means a reaction is active and 0 means it's inactive.
%        cardinality : the max cardinality of iMCS
%        solver: name of the solver either 'gurobi' or 'cplex'
%        tol: tolerance    
%        varargin : objective reaction name(s). If empty biomass is the objective reaction.
%OUTPUT : iMCS : a matrix where one column represents one iMCS. A value of 1 means that the reaction is active 
%         iMCS_card: a vector of size cardinality where iMCS_card(j) is the number of iMCS of cardinality j

clear tic;
clear toc;
tic

%[network, blocked_rxns, change, blocked_mets] = correctModel(network,tol,solver);
[obj_rxns]=get_index_of_rxns(network,varargin{:})

[iMCS, iMCS_card]=hitting_sets(MMB,obj_rxns,cardinality,tol,solver);

%save the iMCS in a structure 
iMCS_res=struct('iMCS',iMCS,'card',iMCS_card);

save('iMCS_res.mat','iMCS_res', '-v7.3');

%check=iMCS_check(IMCS, obj_rxns, new_network,tol,solver);


disp('iMCSs can be found in iMCS_res.mat')

time = toc;
disp(['time to compute iMCS: ', num2str(time)]);

end

