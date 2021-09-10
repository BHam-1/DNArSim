#---------------------------------------------------------------------------#
#	@Authors: Belaid Hamoum, Elsa Dupraz				    #
#	Project : DnarXiv, funded by CominLabs				    #
#	email   : belaid.hamoum@gmail.com				    #
#---------------------------------------------------------------------------#
# 			DNA Data Storage SIMULATOR			    #
#      		     (From Synthesis to Basecalling)			    #
#									    #
#---------------------------------------------------------------------------#

using DelimitedFiles
include("functions.jl")
include("channel.jl")

#parameters
k=parse(Int8,ARGS[4]) #INT [1:10]
nbrSim=parse(Int64,ARGS[1])#INT>0
seq=(readdlm(ARGS[2])) #String: one sequence (with or wothout a header part starting with '>')


include("loadProb.jl") #load the error profile related to k


#simulated reads (channel output)
simSeq=channel(nbrSim,k,seq)



#write simulated reads on a fastq file with no metadata (scores)
outPath=open(ARGS[3],"w")
for iSeq=1:nbrSim
	len=length(simSeq[iSeq])
	write(outPath,">simulation $iSeq \n")
	for j=1:len
		write(outPath,simSeq[iSeq][j])
	end
	write(outPath,"\n+\n")
	write(outPath,"@\n")
end

close(outPath)
