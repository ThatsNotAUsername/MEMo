# -n is mandatory because this is the path of the network. -m is the method choosen it can be elimination (elimination of reversible reactions) or all (we keep all reactions in the network) which is the default method. -s is the solver it can be either cplex (default) or gurobi. -d is the number of decimals we use to solve the problem.
polco_path=/storage/mi/annikaroehl/Tools/polco;
# Absolute path to this script. /home/user/bin/foo.sh
script_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )";
#default parameters. 
solver=cplex
decimal=9

while getopts "n:m:s:d:" opt; do 
 case "$opt" in
    n) network_path=$OPTARG;;
    m) method=$OPTARG;
	    case $method in
		  all) echo "MMBs are going to be computed on the whole cone";;				
		  elimination) echo "After eliminating reversible reactions MMBs will be computed on a reduced cone.";;
		  *) echo "$method does not exist. Please choose elimination/all or don't use this option" ; exit;;
	    esac 
           ;;
    s) solver=$OPTARG;;
    d) decimal=$OPTARG;;			
    *) echo "This option does not exist."; exit;;
 esac
done

# if no network path is given then it fails
if [ -z $network_path ]
then echo 'Failure : There was no network path given' 
     exit;
fi

#This function read a .xml or .mat file containing a metabolic network then correct it and write the text files (matrix.eq and matrix.iq) needed for computing MMBs. We also generate the file that specify what type of reactions we have.
/import/matlabR2015a/bin/matlab -nodisplay -nodesktop -r "pre_processing_mmb('$network_path',1e-$decimal,'$solver'); exit;"
#/import/matlabR2016a/bin/matlab -nodisplay -nodesktop -r "pre_processing_mmb('$network_path',1e-$decimal,'$solver'); exit;"

#If "method" is equal to "elimination" we eliminate reversible reactions with Sage. 1 file is created matrix.eq representing the equalities after the elimination procedure.
if [ "$method" = "elimination" ] 
then time sage matroid_contraction.sage matrix.eq type_rxns.txt $decimal
fi


#compute a generating set of the whole or reduced cone (it depends on the method).
if [ -s matrix.eq ]
then
     #cd /storage/mi/annikaroehl/Tools/polco
     cd $polco_path
     java -jar polco.jar -kind text -eq $script_path/matrix.eq -iq $script_path/matrix.iq -out matlab $script_path/MF_polco.mat
else
     #cd /storage/mi/annikaroehl/Tools/polco
     cd $polco_path
     java -jar polco.jar -kind text -iq $script_path/matrix.iq -out matlab $script_path/MF_polco.mat
fi

#cd /storage/mi/annikaroehl/git-repos/MEMo/MMB\ github
cd $script_path
#Compute the MMB with the generators 
/import/matlabR2016a/bin/matlab -nodisplay -nodesktop -r "from_generators_to_MMB('MF_polco.mat','type_rxns.txt',1e-$decimal); exit;"


#remove the files we created to compute the MMBs
#rm matrix.iq
#rm matrix.eq
#rm type_rxns.txt


