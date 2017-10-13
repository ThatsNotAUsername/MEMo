function [  ] = pre_processing_mmb(network_path,tol,solver)
%This function read a .xml or .mat file containing a metabolic network then correct it and write the text files needed for
%computing MMBs.

%INPUT  network_path : the path of the metabolic network
%       tol : the tolerance used by the solver
%       solver : either 'cplex' or 'gurobi'

if strfind(network_path,'.xml') %the network is a .xml file
    network=readCbModel(network_path);
else net=load(network_path);% the network is .mat file
     names=fieldnames(net);
     network=net.(names{1});
end    
[new_network, blocked_rxns, change, blocked_mets] = correctModel(network,tol,solver);

%save the new_network
save('network.mat','new_network');

n_blocked=sum(blocked_rxns);
n_rev=sum(new_network.rev);
n_irrev=length(new_network.rev)-n_rev;
n_mets=length(new_network.mets);

fprintf('The number of blocked reactions is %d \n',n_blocked)
fprintf('The number of reversible reactions is %d \n',n_rev)
fprintf('The number of irreversible reactions is %d \n',n_irrev)
fprintf('The number of metabolites is %d \n',n_mets)


if not(n_irrev) %there are no irreversible reactions 
    disp('')
    disp('#############################################################')
    disp('')
    error('There are no irreversible reactions. THE NUMBER OF MMBs IS 0')
end

%write three text files : two needed by polco : the stoichiometric matrix (matrix.eq) and the Eirr (Eirr*v>0) matrix (matrix.iq) s.t vi>0 for all reactions i irreversible.
%Also the file type_rxns.txt specify the type of reactions we have in the network 0 :irrev , 1:rev.

%decimal is the number of digits after the decimal point 
decimal=0;
while tol < 1
    decimal=decimal+1;
    tol=10*tol;
end

write_txt_file(new_network,decimal);

end

