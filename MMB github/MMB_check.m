function [ check ] = MMB_check(generators,network,tol)
%This function check if a generator is a MMB. We assume we have generators that only contain irreversible reactions. Then for each generator we remove one active reaction then we solve a LP and see if active reactions have 0 fluxes.  
% Input : generators : each row represent a generator. A value of 1 means
% that the reaction is active 
%          network : a struct with the field .S, .lg , .ub, .rev

clear tic
clear toc
tic

[n_mets n_rxns]=size(network.S);

[n_gen n_rxns]=size(generators);

%vector that specify if a generator is an MMB. 1 means yes 0 means nocheck
check=zeros(1,n_gen);

Model.A=sparse(double(network.S));
Model.lb = network.lb;
Model.ub = network.ub;
vtype=ones(1,n_rxns);
vtype(:)='C';

rev=network.rev;

tosolve = Cplex;
tosolve.Model = Model;
tosolve.Model.obj= zeros(n_rxns,1);
tosolve.Model.lhs = zeros(n_mets,1);
tosolve.Model.rhs = zeros(n_mets,1);
tosolve.Model.ctype = char(vtype);
tosolve.Model.sense='maximize';
    
tosolve.Param.simplex.tolerances.optimality.Cur = tol; % default: 1e-06
tosolve.Param.simplex.tolerances.feasibility.Cur = tol; % default: 1e-06
%tosolve.Param.mip.tolerances.integrality.Cur = tol; % default: 1e-06
tosolve.Param.mip.display.Cur = 'off'; % Solver should not show its iterations
%tosolve.writeModel('cplex_Each_Ind.lp');
tosolve.Param.emphasis.numerical.Cur = true;
    

initial_A=tosolve.Model.A;
initial_rhs=tosolve.Model.rhs;
initial_lhs=tosolve.Model.lhs;
initial_obj=tosolve.Model.obj;  

for i=1:n_gen
    
    gen=generators(i,:);
    
    active_rxns=find(gen);
    
    n_active_rxns=length(active_rxns);
    
    gen_check=zeros(1,n_active_rxns);
    
 % we maximmize the sum of the active rxn fluxes. It has to be 0 to be a MMB
    tosolve.Model.obj(active_rxns)=ones(n_active_rxns,1);

    for r=1:n_active_rxns
                
        irrev=find(~network.rev);
        
        irrev_rxns_not_in_active_rxns=setdiff(irrev,active_rxns);
        
        % add constraints 
    for j=1:n_rxns
        
     
        if rev(j) %the reaction is irreversible
        else vect=zeros(1,n_rxns); %reaction is irreversible
             vect(j)=1;
             tosolve.Model.A(end+1,:)=vect;
             tosolve.Model.lhs(end+1,:)= 0;
             if j==active_rxns(r) | ismember(j,irrev_rxns_not_in_active_rxns)
                tosolve.Model.rhs(end+1,:)=0;
                else tosolve.Model.rhs(end+1,:)=inf;   
             end
        end
    end
    
    %solve a LP
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
            gen_check(r)=1; 
      else
         res=result.objval;
         
         
        if  res<n_active_rxns*tol % all reactions have a 0 flux so the generator is a candidate to be MMB.
            gen_check(r)=1;
        else fprintf('the set %d is not a MMB. When removing rxn %d the objective value is not 0: %f \n',i,active_rxns(r),res)
            find(result.x>tol)
            transpose(result.x)
            break;
        end 
      end
      
      % re-initialize the model for the next reaction r
    tosolve.Model.A=initial_A;
    tosolve.Model.rhs=initial_rhs;
    tosolve.Model.lhs=initial_lhs; 
      
    end
    %re-initialize the objective
    tosolve.Model.obj=initial_obj;

    
    %check if a generator is indeed a mmb
    if gen_check %we have a MMB
       check(i)=1;
    end
end

n_MMB=nnz(check);
if n_MMB==n_gen
    disp('')
    disp('######################')
    disp('All sets are MMBs')
    disp('######################')
    disp('')
else disp('')
    disp('#############################')
    fprintf(' %d sets are not MMBs \n',n_gen-n_MMB )
    disp('#############################')
    disp('')
end


time = toc;
disp(['time to check MMBs: ', num2str(time)]);
    
end

