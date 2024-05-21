# 3D Comparative Assembler using Reference To Assist (3D-CARTA) pipeline

3D Comparative Assembler using Reference To Assist, or 3D-CARTA, is a pipeline can generate a chromosome-length genome assembly fasta from a set of reads, DNA-Seq and/or Hi-C, from an organism by mapping them to an assembly (fragmentary or chromosome-length) of a closely related organism. Similarly to other tools of this type (see, e.g. AMOScmp (Pop et al. 2004) that has served as a primary inspiration for the tool) 3D-CARTA attempts to substitute the traditional overlap-layout-consensus approach to assembly with alignment-layout-consensus. As part of the procedure, the read alignments to the assisting reference (provided in the form of .sam files, typically produced as part of the Juicer workflow) are examined in a base-by-base fashion to detect small-scale differences (single nucleotide polymorphisms, as well as deletions and insertions that are isolated and small as compared to the read length) between the two species. The identified differences are then incorporated into the assisting reference to create a consensus approximation of the organism represented by the read set.

Rearrangements between the two genomes pose the most difficult challenge to comparative assembly using the alignment-layout-consensus strategy. For example, the analysis of read alignment signatures associated with putative rearrangements constitutes the most complex part of the AMOScmp pipeline (Pop et al. 2004). Unlike other comparative assembly tools, 3D-CARTA does not attempt to extract the information necessary to resolve rearrangements from the alignment data. Instead, it relies on Hi-C for this purpose. In particular, it relies on the 3D-DNA/JBAT suite, whose results it takes as input in the form of a .assembly file.

### Overview of the pipeline

`Version date: March 26, 2018`

The pipeline includes the following steps:
	(0) In the preliminary step, the workspace is organized and the assisting reference is preprocessed to create a reference .map file. A .map file is a text file format in which, similar to .mpileup, each line represents a single genomic position.
	(1) During the next step (-S sort) the SAM file(s) are merge-sorted into a master sorted SAM file. The master SAM is then split into manageable chunks. The split size is adjusted dynamically to result in a predefined number of files determined by jobcount parameter (unprompted).
	(3) In the mapping step (-S map) a .map file is constructed from split alignment data to record bases observed at each genomic position. Indels and read break positions with respect to the reference are also recorded.
	(4) The results from individual chunks are merged during the merge step (-S merge). Reference map data is also added during this step. (Including the assisting reference when generating the consensus enhances the accuracy of the assembly in regions of low coverage, or regions covered only by a single read, helping distinguish the true base from contaminants and sequencing errors.) Optionally, read clips that originate from restriction sequence-based cutting in the Hi-C protocol can be filtered during this step.
	(5) The next step (-S reduce) generates a preliminary consensus from the map.
	(6) The pipeline includes an optional step to reexamine indel and clip features that were not incorporated during preliminary data consolidation and reassemble the corresponding regions (-S reassemble).
	(7) In the last step (-S finalize) the rearrangement file (provided via -a|--assembly flag) is incorporated and the final FASTA is constructed that reflects both the local consensus and the rearrangement data.

#### This software is distributed under The MIT License (MIT).

#### Prerequisites
- `Bash >=4`
- `GNU Awk >=4.0.2`
- `GNU coreutils sort >=8.11`
- `GNU Parallel >=20150322` – highly recommended to increase performance

#### Installation
There is no need for installation.

#### USAGE
```sh
--------------------
Usage: ./run-3d-carta.sh [options] <path_to_reference_fasta> <path_to_sam_file(s)>

ARGUMENTS:

path_to_reference				Path to an assisting species reference fasta.
path_to_sam_file(s)        		Path to sam file(s) containing alignments of data generated using a targe species sample to an assisting species reference.

MAIN OPTIONS:

-h|--help                      	
								Shows this help.
-a|--assembly		path_to_assembly_file	
								Path to assembly file (typically from 3D-DNA and JBAT), describing a set of large-scale assembly actions (cut, anchor, order and orient) that are required to create a genome for species of interest from the assisting reference.
-c|--consensus		consensus_threshold 
								The number of bases in the alignment data accepted to call a basepair in the final sequence, i.e. how many times the base needs to be seen across the read set (and reference if c>1) to constitute a valid call. Default: 1.
-q|--mapq			mapq_threshold  
								Minimum mapping quality of the read to be considered in the recontruction process. Default: 0.
-r|--restriction	restriciton_sequence
								Restriction sequence for the enzyme. Expected if using Hi-C read alignments as input. Default: GATC.
								

SUPPLEMENTARY OPTIONS:

** worflow **

-t|--threads		number_of_threads
								Specify the number of threads to run.
-s|--stage			stage
								Start from a particular stage. Can be sort, map, merge, reduce, reassemble and finalize.
--single-stage
								Run only one stage, whichever is passed on with -s|--stage (if -s|--stage is not triggered only preliminary stage is run).

** tuning **

--trust-reference				
								Consider reference as contributing to consensus [applicable only for c>1, if c=1 trust_reference is ignored].


** tuning **

--trust-reference                               
                                                                Consider reference as contributing to consensus [applicable only for c>1, if c=1 trust_reference is ignored].
```sh
