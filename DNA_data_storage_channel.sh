#---------------------------------------------------------------------------#
#	@Authors: Belaid Hamoum, Elsa Dupraz				    #
#	Project : DnarXiv, funded by CominLabs				    #
#	email   : belaid.hamoum@gmail.com				    #
#---------------------------------------------------------------------------#
# 			DNA Data Storage SIMULATOR			    #
#      		     (From Synthesis to Basecalling)			    #
#									    #
#---------------------------------------------------------------------------#

	####
	#->Parameters:
	#	-$1:Seq path  to simulate
	#	-$2:Nbr of reads to simulate
	#	-$3:Simulated sequences output path
	#	-$4:Channel memory length, k=6 by default [RECOMENDED] 
	###

cd $(dirname $0)
prevDir=$(pwd)


	inPath=""
	k=6
	nbrRead=0
	outPath=""

	while [ "$1" != "" ]; do
	    case $1 in
		-i | --inputSeq )
		    shift
		    inPath=$1
		;;
		-o | --outPath )
		    shift
		    outPath=$1  
		;;            
		-k | --memLen )    
		    shift
		    k=$1
		;;
		-n | --nbrRead )   
		    shift
		    nbrRead=$1
		;;

	    esac
	    shift
	done

	#inputs checking

	if [ $k == "" ]
	then
		k=6
	
	else
		int='^[1-9][0-9]*$'
		if ! [[ $k =~ $int ]]
		then
		   echo "error: -k [1:10] (pick up a value between 1 and 10), [k=6 is recommended]"
		   exit 1
		elif [[ $k -gt "10" ]]
		then
		   echo "error: -k [1:10] (pick up a value between 1 and 10), [k=6 is recommended]"
		   exit 1
		fi
	fi

	

	int='^[0-9]+$'
	if ! [[ $nbrRead =~ $int ]] || [ $nbrRead -eq "0" ]
	then
		echo "error: value of 'n' should be greater than 0 -n [1:N]"
		exit 1
	
	fi
	


	if [ -f "$prevDir/$inPath" ]
	then
		inPath=$prevDir/$inPath
	elif [ -f "$inPath" ]
	then
		inPath=$inPath
	else
		echo "error: input path file  is incorrect"
		exit 1
	fi

	
	


	if ! touch $prevDir/$outPath > /dev/null 2>&1 
	then
		if ! touch $outPath > /dev/null 2>&1 
		then
			echo "error: incorrect output path"
			exit 1
		else
			outPath=$outPath
		fi

	else
		outPath=$prevDir/$outPath
	fi




	#launch
	
	echo -n "-> Launching... "

	cd ./simulator

		julia simulator_v1.jl $nbrRead $inPath $outPath $k

	cd ..


cd $prevDir
