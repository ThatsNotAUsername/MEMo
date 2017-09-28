function [ ] = write_txt_file(network,decimal)
%write three text files : two needed by polco : the stoichiometric matrix
%(matrix.eq) and the Eirr matrix (matrix.iq) s.t Eirr*v>0.
%Also the file type_rxns.txt will contain the type of reactions we have in
%the network 0 :irrev , 1:rev and 2:blckd.

%INPUT : network :network struct with the field S and rev.
%        decimal : number of digits after the decimal point used to save the stoichiometric
%        matrix in the text file matrix.eq
        

S=network.S;
stoichiometricMatrix=full(S);

[n_mets,n_rxns]=size(S);


fid1 = fopen('/home/mi/triou/Documents/matrix.eq','w');

%copy the S matrix in a text file called matrix.eq
for i=1:n_mets
    for j=1:n_rxns
    fprintf(fid1,['%.' num2str(decimal) 'f %.' num2str(decimal) 'f'],stoichiometricMatrix(i,j));
    end
    fprintf(fid1,'\n');
end 

fclose(fid1);

fid2=fopen('/home/mi/triou/Documents/matrix.iq','w');
fid3=fopen('/home/mi/triou/Documents/type_rxns.txt','w');



for i=1:n_rxns
    if network.rev(i) % the reaction is reversible there is no constraint
        fprintf(fid3,'%d %d',1);
    else  fprintf(fid3, '%d %d',0);% the reaction is irreversible
        for j=1:(i-1)
            fprintf(fid2,'%d %d',0); 
        end    
        fprintf(fid2,'%d %d',1);
        for j=(i+1):n_rxns
            fprintf(fid2,'%d %d',0); 
        end  
     fprintf(fid2,'\n',0); 

    end
end

fclose(fid2);
fclose(fid3);

end

