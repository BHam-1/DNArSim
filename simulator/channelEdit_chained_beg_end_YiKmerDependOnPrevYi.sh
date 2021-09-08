#---------------------------------------------------------------------------#
#	@Authors: Belaid Hamoum, Elsa Dupraz				    #
#	Project : DnarXiv, funded by CominLabs				    #
#	email   : belaid.hamoum@gmail.com				    #
#---------------------------------------------------------------------------#
# 			DNA Data Storage SIMULATOR			    #
#      		     (From Synthesis to Basecalling)			    #
#									    #
#---------------------------------------------------------------------------#


#$1:reference sequence, $2: number of sequences to generate (simulate), $3:output folder, $4:memory length



	
	refFile=tmpRefFile.txt
	
	echo $1 > $refFile

	n=$2 #nbr of seq to simulate


	k=6

	if [ ! -z "$4" ]
	then
		k=$4

	fi

	simOut=$3

	echo $n
	echo $refFile
	echo $simOut
	echo $k

	julia simulator_v1.jl $n $refFile $simOut $k


cd $prevDir
