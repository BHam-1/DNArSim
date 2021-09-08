# DNA_Data_Storage_Channel_Model

Simulator to simulate the entire DNA Data Storage process (from synthesis to basecalling)


This is the first version of the Channel Model with Memory for DNA Data Storage with Nanopore Sequencing based on Julia (https://julialang.org/downloads/).

After setting Julia correctly, simulations can be launched using next command:

> ./DNA_data_storage_channel.sh -i ./example/ref.txt -n 100 -o ./example/sim.fastq -k 6
	
**Parameters**:

 * **-i**:  Path to the input sequence (should be on fasta format) to simulate. An example of such sequence is available on "example" folder with the name 'ref.txt'. [required]
 * **-n**:  Number of reads to simulate. [required]
 * **-o**:  Simulated sequences output path. Will be presented in a fastq format (without included scores). [required]
 * **-k**:  Channel memory length.  fixed to k=6 by default [recommended]
	
More documentation will be added as soon as possible. Meanwhile feel free to reach me by e-mail (**belaid.hamoum@gmail.com**) for more details.
