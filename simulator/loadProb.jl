#---------------------------------------------------------------------------#
#	@Authors: Belaid Hamoum, Elsa Dupraz				    #
#	Project : DnarXiv, funded by CominLabs				    #
#	email   : belaid.hamoum@gmail.com				    #
#---------------------------------------------------------------------------#
# 			DNA Data Storage SIMULATOR			    #
#      		     (From Synthesis to Basecalling)			    #
#									    #
#---------------------------------------------------------------------------#


function readProbFile(path)

	return Float32.(readdlm(path)[1])

end

root="./probEdit/k$(k)"

BegProbEdit=readProbFile("$(root)/BegErrByPosAvg.txt")
BegProbIns=readProbFile("$(root)/BegInsByPosAvg.txt")
BegProbDel=readProbFile("$(root)/BegDelByPosAvg.txt")
BegProbSubst=readProbFile("$(root)/BegMisAvg.txt")


MidProbEditA=readProbFile("$(root)/AErrAvg.txt")
MidProbInsA=readProbFile("$(root)/AInsAvg.txt")
MidProbDelA=readProbFile("$(root)/ADelAvg.txt")
MidProbSubstA=readProbFile("$(root)/AMisAvg.txt")


MidProbEditC=readProbFile("$(root)/CErrAvg.txt")
MidProbInsC=readProbFile("$(root)/CInsAvg.txt")
MidProbDelC=readProbFile("$(root)/CDelAvg.txt")
MidProbSubstC=readProbFile("$(root)/CMisAvg.txt")


MidProbEditG=readProbFile("$(root)/GErrAvg.txt")
MidProbInsG=readProbFile("$(root)/GInsAvg.txt")
MidProbDelG=readProbFile("$(root)/GDelAvg.txt")
MidProbSubstG=readProbFile("$(root)/GMisAvg.txt")


MidProbEditT=readProbFile("$(root)/TErrAvg.txt")
MidProbInsT=readProbFile("$(root)/TInsAvg.txt")
MidProbDelT=readProbFile("$(root)/TDelAvg.txt")
MidProbSubstT=readProbFile("$(root)/TMisAvg.txt")


EndProbEdit=readProbFile("$(root)/EndErrByPosAvg.txt")
EndProbIns=readProbFile("$(root)/EndInsByPosAvg.txt")
EndProbDel=readProbFile("$(root)/EndDelByPosAvg.txt")
EndProbSubst=readProbFile("$(root)/EndMisAvg.txt")


subList=["A2C" "A2G" "A2T";
 "C2A" "C2G" "C2T";
 "G2A" "G2C" "G2T";
 "T2A" "T2C" "T2G"]


Yi=["M" "D" "I" "S"]



transProb=[]
for i=1:length(subList[:,1])
	tmpProb=zeros(Float32,0)
	for j=1:length(subList[1,:])
		avgTmp=readProbFile("$(root)/$(subList[i,j])_Avg.txt");
		push!(tmpProb,avgTmp);
	end
	push!(transProb,tmpProb)
end

#I should get thes RESULTS DIRECTLY DURING the training !!!!!----- RESOLVE THIS!!!!!!!!!!!!
for i=1:length(subList[:,1])
	transProb[i]=transProb[i]/sum(transProb[i])
end

rangeTransProb=[]

for i=1:length(subList[:,1])
	tmpRangeLine=[]
	for j=1:length(subList[1,:])
		if(j>1)
			push!(tmpRangeLine,tmpRangeLine[j-1]+transProb[i][j]);
		else
			push!(tmpRangeLine,transProb[i][j]);
		end
	end
	push!(rangeTransProb,tmpRangeLine)
end



EditYiKmerprevYi=[]
mapKmerPrevYi=[]

for e1=1:4
	if(isfile("$(root)/KmerYi_prevYi$(Yi[e1])_RatesAvg.txt") )
		#kmer Ins Del Subst Err Match
		push!(EditYiKmerprevYi,readdlm("$(root)/KmerYi_prevYi$(Yi[e1])_RatesAvg.txt"));
		nKmer=length(EditYiKmerprevYi[e1][:,1]) #nbr of entries


		push!(mapKmerPrevYi,Dict())

		for i=1:nKmer
			mapKmerPrevYi[e1][EditYiKmerprevYi[e1][i,1]]=i
		end

	else
		push!(mapKmerPrevYi,Dict())
		push!(EditYiKmerprevYi,[]);
	end

end


#index: 1-> prev Match, 2-> prev Del, 3-> prev Ins, 4-> prev Subst,
kmerInsLenProbPrevYi=[]
mapKmerInsLenPrevYi=[]
#Yi-1= [Match[Yi=I,D,S,E,M],Del[Yi=I,D,S,E,M],Ins[Yi=I,D,S,E,M],Sub[Yi=I,D,S,E,M]]
for e1=1:4
	if(isfile("$(root)/KmerInsLen_prevYi$(Yi[e1])_RatesAvg2.txt") )
		kmerInsLenProbTmp=readdlm("$(root)/KmerInsLen_prevYi$(Yi[e1])_RatesAvg2.txt");



		push!(mapKmerInsLenPrevYi,Dict())
		push!(kmerInsLenProbPrevYi,[])
		nKmer=length(kmerInsLenProbTmp[:,1])

		for i=1:nKmer
			mapKmerInsLenPrevYi[e1][kmerInsLenProbTmp[i,1]]=i
			tmpLine=parse.(Float32,split(join(kmerInsLenProbTmp[i,2:end]),","))

			for j=1:length(tmpLine)
				if(j>1)
						tmpLine[j]+=tmpLine[j-1]
				end
			end
			push!(kmerInsLenProbPrevYi[e1],tmpLine)
		end
		kmerInsLenProbTmp=Nothing;



	else
		push!(mapKmerInsLenPrevYi,Dict())
		push!(kmerInsLenProbPrevYi,[])
	end


end





#probDel of current k-mer knowing that a delete occurs previously (Yi(kmeri) depend on Y-1)
DelProb_YiprevDel=readdlm("$(root)/KmerDel_delPrevRates.txt");
mapDelProb_YiprevDel=Dict()
nKmer=length(DelProb_YiprevDel[:,1]) #nbr of entries
for i=1:nKmer
	mapDelProb_YiprevDel[DelProb_YiprevDel[i,1]]=DelProb_YiprevDel[i,2]
end


probInsLenBeg=Float32.(readdlm("$(root)/insLenBegRates.txt")[1,1:end-1])
probInsLenMidA=Float32.(readdlm("$(root)/AinsLenMidRates.txt")[1,1:end-1])
probInsLenMidC=Float32.(readdlm("$(root)/CinsLenMidRates.txt")[1,1:end-1])
probInsLenMidG=Float32.(readdlm("$(root)/GinsLenMidRates.txt")[1,1:end-1])
probInsLenMidT=Float32.(readdlm("$(root)/TinsLenMidRates.txt")[1,1:end-1])
probInsLenEnd=Float32.(readdlm("$(root)/insLenEndRates.txt")[1,1:end-1])

sumBeg=BegProbIns+BegProbDel+BegProbSubst
rangeInsMaxBeg=BegProbIns/sumBeg;
rangeDelMaxBeg=rangeInsMaxBeg+(BegProbDel/sumBeg);
rangeSubstMaxBeg=rangeDelMaxBeg+(BegProbSubst/sumBeg);

sumMidA=MidProbDelA+MidProbSubstA
rangeDelMaxMidA=(MidProbDelA/sumMidA);
rangeSubstMaxMidA=rangeDelMaxMidA+(MidProbSubstA/sumMidA);

sumMidC=MidProbDelC+MidProbSubstC
rangeDelMaxMidC=(MidProbDelC/sumMidC);
rangeSubstMaxMidC=rangeDelMaxMidC+(MidProbSubstC/sumMidC);

sumMidG=MidProbDelG+MidProbSubstG
rangeDelMaxMidG=(MidProbDelG/sumMidG);
rangeSubstMaxMidG=rangeDelMaxMidG+(MidProbSubstG/sumMidG);

sumMidT=MidProbDelT+MidProbSubstT
rangeDelMaxMidT=(MidProbDelT/sumMidT);
rangeSubstMaxMidT=rangeDelMaxMidT+(MidProbSubstT/sumMidT);

sumEnd=EndProbIns+EndProbDel+EndProbSubst
rangeInsMaxEnd=EndProbIns/sumEnd;
rangeDelMaxEnd=rangeInsMaxEnd+(EndProbDel/sumEnd);
rangeSubstMaxEnd=rangeDelMaxEnd+(EndProbSubst/sumEnd);

lenInsBeg=length(probInsLenBeg);
rangeLenInsBeg=[]
for j=1:lenInsBeg
	if(j>1)
		push!(rangeLenInsBeg,rangeLenInsBeg[j-1]+probInsLenBeg[j])
	else
		push!(rangeLenInsBeg,probInsLenBeg[1])
	end
end


lenInsMidA=length(probInsLenMidA);
rangeLenInsMidA=[]
for j=1:lenInsMidA
	if(j>1)
		push!(rangeLenInsMidA,rangeLenInsMidA[j-1]+probInsLenMidA[j])
	else
		push!(rangeLenInsMidA,probInsLenMidA[1])
	end
end

lenInsMidC=length(probInsLenMidC);
rangeLenInsMidC=[]
for j=1:lenInsMidC
	if(j>1)
		push!(rangeLenInsMidC,rangeLenInsMidC[j-1]+probInsLenMidC[j])
	else
		push!(rangeLenInsMidC,probInsLenMidC[1])
	end
end

lenInsMidG=length(probInsLenMidG);
rangeLenInsMidG=[]
for j=1:lenInsMidG
	if(j>1)
		push!(rangeLenInsMidG,rangeLenInsMidG[j-1]+probInsLenMidG[j])
	else
		push!(rangeLenInsMidG,probInsLenMidG[1])
	end
end

lenInsMidT=length(probInsLenMidT);
rangeLenInsMidT=[]
for j=1:lenInsMidT
	if(j>1)
		push!(rangeLenInsMidT,rangeLenInsMidT[j-1]+probInsLenMidT[j])
	else
		push!(rangeLenInsMidT,probInsLenMidT[1])
	end
end

lenInsEnd=length(probInsLenEnd);
rangeLenInsEnd=[]
for j=1:lenInsEnd
	if(j>1)
		push!(rangeLenInsEnd,rangeLenInsEnd[j-1]+probInsLenEnd[j])
	else
		push!(rangeLenInsEnd,probInsLenEnd[1])
	end
end
