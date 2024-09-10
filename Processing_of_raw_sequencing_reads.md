#Processing of raw sequencing reads

After obtaining raw sequencing reads, they were mapped with a combined reference genome of human GRCh38 and *E. intestinalis* ATCC 50506 using Cell Ranger `count` that was installed on HPC. 

##**Procedure**

Log in to BigPurple 

`ssh -Y jaroep01@bigpurple.nyumc.org`

Connect to cpu node before calling Cell Ranger
 
`srun -p cpu_dev -n 3 --mem-per-cpu=20G --time=2:00:00 --x11 --pty bash`

Load Cell Ranger 4.0.0 and bcl2fastq 2.2.0

`module load cellranger/4.0.0`

`module load bcl2fastq`

To demultiplex the dataset containing hashing oligos (HTO). We need to create 1) library.csv and 2) feature_reference.csv that looked like these

![1](https://github.com/pjaroenlak/scRNAseq_E_intestinalis/blob/main/Images/Screen_Shot_2020-11-09_at_6.25.08_PM.png)

![2](https://github.com/pjaroenlak/scRNAseq_E_intestinalis/blob/main/Images/Screen_Shot_2020-11-09_at_6.38.38_PM.png)

Next, we are going to make the combined reference genome. We will run this using SBATCH to make it runs eventhough we lost internet connection

First, create the text file the commands we are going to run

`touch queue2.txt`

`vi queue2.txt`


	#!/bin/bash

	#SBATCH --partition=gpu4_medium
	#SBATCH --job-name=Michael_job
	#SBATCH --ntasks=3
	#SBATCH --cpus-per-task=2
	#SBATCH --output=log_%A_%a.out
	#SBATCH --error=log%A_%a.err
	#SBATCH --mem=250G
	#SBATCH --time=48:00:00

	module load cellranger/4.0.0
	module load bcl2fastq

	cellranger mkref --genome=Homo_Ei_ref_v3 --fasta=./GRCh38-2020-A_build/	Human_Eintestinalis_Genome_v3.fa --genes=./GRCh38-2020-A_build/	Human_Eintestinalis_v3.filtered.gtf


`sbatch queue2.txt `

Next, we are going to map the fastq reads to the combined reference genomes. 

`touch queue.txt`
`vi queue.txt`

	#!/bin/bash

	#SBATCH --partition=gpu4_medium
	#SBATCH --job-name=Michael1
	#SBATCH --ntasks=3
	#SBATCH --cpus-per-task=2
	#SBATCH --output=log_%A_%a.out
	#SBATCH --error=log%A_%a.err
	#SBATCH --mem=250G
	#SBATCH --time=48:00:00

	module load cellranger/4.0.0
	module load bcl2fastq

	cellranger count --id Macrophages_Ei_v3 --libraries /gpfs/data/bhabhaekiertlabs/home/jaroep01/MichaelPJ/scRNA-seq/mapping_v3/library_macrophage_Ei.csv --transcriptome /gpfs/data/bhabhaekiertlabs/home/jaroep01/MichaelPJ/scRNA-seq/reference-genome/Homo_Ei_ref_v3 --feature-ref /gpfs/data/bhabhaekiertlabs/home/jaroep01/MichaelPJ/scRNA-seq/mapping_v3/feature_ref_macrophage_Ei.csv

`run sbatch`

`sbatch queue.txt`

Output files are stored at 

`/gpfs/data/bhabhaekiertlabs/home/jaroep01/MichaelPJ/scRNA-seq/mapping_v3/Macrophages_Ei_v3/outs`