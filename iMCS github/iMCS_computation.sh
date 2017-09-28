# -n is mandatory because this is the path of the network. -m is the method choosen it can be elimination (elimination of reversible reactions) or all (we keep all reactions in the network) which is the default method. -c is the maximum cardinality of iMCSs we are looking for (default is 1e09). -s is the solver it can be either cplex (default) or gurobi. -d is the number of decimal we use solve the problem. -r the objective reaction default is BIOMASS

#default parameters. 
solver=cplex
decimal=9
cardinality=100
objective=BIOMASS

while getopts "n:m:c:s:d:r:" opt; do 
 case "$opt" in
    n) network_path=$OPTARG;;
    m) method=$OPTARG;
	    case $method in
		  all) echo "MMBs are going to be computed on the whole cone";;				
		  elimination) echo "After eliminating reversible reactions MMBs will be computed on a reduced cone.";;
		  *) echo "$method method does not exist. Please choose elimination/all or  don't use this option" ; exit;;
	    esac 
           ;;
    c) cardinality=$OPTARG;;
    s) solver=$OPTARG;
		case $solver in
		  cplex) echo "Cplex solver was chosen";;				
		  gurobi) echo "Gurobi solver was chosen";;
		  *) echo "$solver solver does not exist. Please choose cplex/gurobi or don't use this option" ; exit;;
	    esac 
           ;;
    d) decimal=$OPTARG;;
    r) objective=$OPTARG;;			
    *) echo "This option does not exist."; exit;;
 esac
done

# if there is no path of network given then it fails
if [ -z $network_path ]
then echo 'Failure : There was no path of network given' 
     exit;
fi
            

#This function read a .xml or .mat file containing a metabolic network then correct it and write the text files (matrix.eq and matrix.iq) needed for computing MMBs. We also generate the file that specify what type of reactions we have.
/import/matlabR2016a/bin/matlab -nodisplay -nodesktop -r "pre_processing_mmb('$network_path',1e-$decimal,'$solver'); exit;"


#If "method" is equal to "elimination" we eliminate reversible reactions with Sage. 1 file is created matrix.eq representing the equalities after the elimination procedure.
if [ "$method" = "elimination" ]  
then sage matroid_contraction.sage matrix.eq type_rxns.txt $decimal
fi



#compute a generating set of the whole or reduced cone (it depends on the method).
if [ -s matrix.eq ]
then
     java -jar polco.jar -kind text -eq matrix.eq -iq matrix.iq -out matlab MF_polco.mat
else
     java -jar polco.jar -kind text -iq matrix.iq -out matlab MF_polco.mat
fi


#Compute the MMB with the generators and then the iMCS. MMB are saved in a file MMB.mat and iMCS in a file iMCS_res.mat
/import/matlabR2016a/bin/matlab -nodisplay -nodesktop -r "MMB=from_generators_to_MMB('MF_polco.mat','type_rxns.txt',1e-$decimal); network=load('network.mat'); network=network.new_network; iMCS_computation( network,MMB,$cardinality,1e-$decimal,'$solver','$objective'); exit;"

#remove the files we created to compute the MMBs
rm matrix.iq
#rm matrix.eq
rm type_rxns.txt
rm MF_polco.mat

