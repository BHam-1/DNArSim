#---------------------------------------------------------------------------#
#	@Authors: Belaid Hamoum, Elsa Dupraz				    #
#	Project : DnarXiv, funded by CominLabs				    #
#	email   : belaid.hamoum@gmail.com				    #
#---------------------------------------------------------------------------#
# 			DNA Data Storage SIMULATOR			    #
#      		     (From Synthesis to Basecalling)			    #
#									    #
#---------------------------------------------------------------------------#



function randomBase()
	base=rand((0,1,2,3))

	if(base==0)
		return "A"
	elseif(base==1)
		return "C"
	elseif(base==2)
		return "G"
	else
		return "T"
	end
end

function insertB(seqTmp,rangeLenIns,lenIns)
	r3=rand()
	len=1;

	#length of the insertion
	for j=1:lenIns

		if (r3<=rangeLenIns[j])
			len=j;
			break;
		end

	end

	for j=1:len #insert after current pos
		base=randomBase() #pickup randomly a base to insert
		push!(seqTmp,base)
	end

end


function deleteB(seqTmp)
	len=1; #length of del is 1 (bursts are taken into account in another way)
	#Do nothing => deletion
end

function substB(seqTmp,rangeTransProb,currBase)
	r=rand()
	base=""
	if(currBase=='A')
		for j=1:length(rangeTransProb[1])
			if (r<=rangeTransProb[1][j])
				base=split(subList[1,j],"2")[2] 
				break;
			end
		end
	elseif(currBase=='C')
		for j=1:length(rangeTransProb[2])
			if (r<=rangeTransProb[2][j])
				base=split(subList[2,j],"2")[2] 
				break;
			end
		end
	elseif(currBase=='G')
		for j=1:length(rangeTransProb[3])
			if (r<=rangeTransProb[3][j])
				base=split(subList[3,j],"2")[2] 
				break;
			end
		end
	elseif(currBase=='T')
		for j=1:length(rangeTransProb[4])
			if (r<=rangeTransProb[4][j])
				base=split(subList[4,j],"2")[2]
				break;
			end
		end
	end
	push!(seqTmp,base)
end


function matchB(seqTmp,currBase)
	push!(seqTmp,currBase)
end

