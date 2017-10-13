import numpy as np
import sys


def eliminate_reversible(S,rev_list,decimal):

#This function writes one file 'matrix.eq' that contains the equalities in the system after elimination of reversible reactions.
#INPUT : S, the stoichiometric matrix
#	 rev_list, a vector s.t rev_list(i)= 1 if a reaction is reversible and 0 if it's irreversible.

	#M is the linear matroid generated by the stoichiometric matrix S
	M=Matroid(matrix=S)

	#Contract the reversible reactions
	M=M/rev_list


	# the groundset_list is not sorted. We concatenate the groundset_list and the representation of the matroid M after contraction. Then we sort the matrix according to the groundset_list.
	mat_to_sort=np.concatenate(([M.groundset_list()],M.representation()),axis=0)

	mat_to_sort=mat_to_sort[:,mat_to_sort[0,:].argsort()]

	sorted_mat=mat_to_sort[1:len(mat_to_sort),:]

	eq=np.zeros([len(sorted_mat),len(data_sto[0])]) 

	#the matrix eq has n_rxns=n_irrev+n_rev column. The reversible columns are 0, the irreversible columns are filed with sorted_mat
	groundset_list=sorted(M.groundset_list())

	for i  in range(0,len(sorted_mat)):
		eq[i,groundset_list]=sorted_mat[i,:]

	#the equalities are stored in .eq file
	np.savetxt('matrix.eq',eq,delimiter=' ',fmt='%.'+str(decimal)+'f');


#read the file that contains the stoichiometric matrix
data_sto=np.genfromtxt(str(sys.argv[1]),delimiter='')
S=Matrix(QQ,data_sto)


#read the file that specify if a reaction is reversible or irreversible
data_rxns=np.genfromtxt(str(sys.argv[2]),delimiter='')
rev_list=np.where(data_rxns==1)[0].tolist()

#read the number of decimals
decimal=sys.argv[3]


#eliminate reversible reactions
eliminate_reversible(S,rev_list,decimal)

