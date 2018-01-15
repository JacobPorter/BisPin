nohup ~/BisPin/BisPin_align.py -N -i 1 -n 2 -x ~/bfast-gap-0.1.0/exponential_scoring.txt  -g -P /home/jsporter/bfast-gap-0.1.0/bfast-gap-original/bfast-gap /research/jsporter/Data/genome/GRCh38.p9/BisFAST/GRCh38.p9.fa ./test.original.sam /research/jsporter/Data/reads/Human/IonTorrent/Simulation/test.2000.fastq &> ./test.original.out &

nohup ~/BisPin/BisPin_align.py -N -i 1 -n 2 -x ~/bfast-gap-0.1.0/exponential_scoring.txt  -g -P /home/jsporter/bfast-gap-0.1.0/bfast-gap-lookup/bfast-gap /research/jsporter/Data/genome/GRCh38.p9/BisFAST/GRCh38.p9.fa ./test.lookup.sam /research/jsporter/Data/reads/Human/IonTorrent/Simulation/test.2000.fastq &> ./test.lookup.out &

nohup ~/BisPin/BisPin_align.py -N -i 1 -n 2 -x ~/bfast-gap-0.1.0/exponential_scoring.txt  -g -P /home/jsporter/bfast-gap-0.1.0/bfast-gap/bfast-gap /research/jsporter/Data/genome/GRCh38.p9/BisFAST/GRCh38.p9.fa ./test.runCounts.sam /research/jsporter/Data/reads/Human/IonTorrent/Simulation/test.2000.fastq &> ./test.runCounts.out &
