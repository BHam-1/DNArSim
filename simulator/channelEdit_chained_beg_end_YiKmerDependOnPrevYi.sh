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



prevDir=$(pwd)

cd $(dirname $0)
	
	simOut=""
	# check if absolute or relative
	if [ -f $prevDir/$3 ] 
	then	
		simOut=$prevDir/$3
	elif [ -f $3 ]
	then
		simOut=$3
	fi
	
	refFile=tmpRefFile.txt
	
	echo $1 > $refFile

	n=$2; #nbr of seq to simulate


	k=6

	if [ ! -z "$4" ]
	then
		k=$4

	fi
	

	julia simulator_v1.jl $n $refFile $simOut $k


cd $prevDir
