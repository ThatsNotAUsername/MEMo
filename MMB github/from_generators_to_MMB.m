function [ MMB ] = from_generators_to_MMB(generators_file,rev_file,tol)
%This function creates one text file containing the MMB
%INPUT : generators_file: path of the mat file containg the generators (field .R contains the generators in column)
%        rev_file: path of the file specifying what type of reaction we have
%        tol : tolerance all numbers less than tol are considered as
%        equal to 0

%OUTPUT : MMB: matrix of the MMB where each row represent an MMB and each
%              column a reaction. A value of 1 in row i column j means that reaction j is
%              active in MMB i. Otherwise we put 0.


%read the generators
polco_output=load(generators_file);

ER=transpose(polco_output.R);


[n_ER n_rxns]=size(ER);

%read the file that specify the status of each reaction 
fid=fopen(rev_file,'r');

type_rxns=fscanf(fid,'%d');

irrev=find(type_rxns==0);

fclose(fid);

%create the file where we will store the MMBs
fid2=fopen('MMB_polco.txt','w');

MMB=[];
vect_mmb=zeros(1,n_rxns);
n_MMB=0;

% a MMB only contains irreversible reactions so we only look at the irreversible reaction with a non zero flux. 

for i=1:n_ER
    active_rxns=find(ER(i,:)>tol);
    mmb=intersect(active_rxns,irrev);
    if not(isempty(mmb))
    vect_mmb(mmb)=ones(1,length(mmb));    
    MMB=[MMB ; vect_mmb];
    n_MMB=n_MMB+1;
    vect_mmb=zeros(1,n_rxns);
  
    %write in a text file the MMBs
%     fprintf(fid2,'%d ',mmb);
%  
%      fprintf(fid2,'\n');
%      fprintf(fid2,'\n');
    end
end

save('MMB.mat','MMB');

fclose(fid2);

fprintf('The number of MMBs is %d \n',n_MMB)
fprintf('MMBs can be found in MMB.mat \n')
   
end

