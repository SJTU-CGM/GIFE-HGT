# GIFE-HGT: a Genome-wide Identifier For Eukaryote Horizontal Gene Transfers
GIFE-HGT is a systematic and efficient method for genome-wide identification of horizontal gene transfers (HGTs) in eukaryotic genomes. GIFE-HGT incorporates an accelerated mode specifically designed for analyzing large genomes, which maintains high accuracy while significantly reducing computational requirements. Unlike conventional approaches that are limited to protein-coding regions, our method enables the discovery of novel HGTs in non-coding and regulatory regions. By implementing a strategic taxonomic classification system consisting of self-group (SG), closely related group (CRG), and distantly related group (DRG), GIFE-HGT effectively differentiates between ancient and recent HGT events.

## Installation
### Requirements

- Perl 5.26.3 or up
- Python 3.11
- Software needed: [lastz][1], [TRF][2], [cd-hit][3], [Bowtie2][4], [BLAST][5],[samtools][6], [MAFFT][7], [trimAl][8], [IQ-TREE][9].
  You can also refer to this website (https://cgm.sjtu.edu.cn/GIFE-HGT/install.html) to install the above software.

### Installation procedures
#### 1. Download the GIFE-HGT from github
```
$ git clone https://github.com/SJTU-CGM/GIFE-HGT.git
```
#### 2. Add bin/ to PATH and add lib/ to LD_LIBRARY_PATH
```
$ export PATH=$PATH:/path/to/GIFE-HGT/bin/:
$ export PERL5LIB=$PERL5LIB:/path/to/GIFE-HGT/lib/:
$ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/GIFE-HGT/lib/:
```
#### 3. Test if GIFEHGT is installed successfully: using `GIFEHGT`. 
If you see the following content, congratulations! GIFE-HGT is successfully installed. If not, see if all the requirements are satisfied.
```
Usage: GIFEHGT <command> ...

Available commands:
	kmerFilter	Select most different fragments with target genome by calculating the distance of kmer frequecies
	splitDB		Split the genomic database into three groups: SG, CRG and DRG
	seqAlign	Sequence aligment using LASTZ
	screenHGT	Screen potential HGTs using sequence alignment results
	WGSValidate	Validate potential HGTs using WGS datasets
	conPhyTree	Construct sequence phylogenetic tree to validate HGTs
```
The usage information for each command can be shown with the --help or -h option after each command name. 

## Quick start
### Example data
We use the identification of HGTs in one chromosome of yeast with some bacteria as the example data.   
Please download from here (http://cgm.sjtu.edu.cn/GIFE-HGT/data/example.tar.gz) and decompress it:
```
$ wget http://cgm.sjtu.edu.cn/GIFE-HGT/data/example.tar.gz
$ tar zxvf example.tar.gz & cd example
```
The dataset includes 
-  `data/fa/sacCer3_NC_001144.5.fa`: sequence of one chromosome of yeast
-  `data/db/refseq/*`: genomic databases including 10 fungi and 13 bacteria
-  `data/db/*txt`: taxonomic information of all organisms in genomic datasets
-  `data/db/WGSdata/*`:  Whole genomic sequencing datasets including several bacteria
-  `data/mitochondria_chloroplast/*`: mitochondrial sequence of yeast
-  `data/rmsk/*`: repeat sequence annotation file of yeast
-  `data/singlecopy/*`: single copy genes common in all eukaryotes.
### Manual (`Example command` part uses example data mentioned above)
#### 1. `kmerFilter`: Select most different fragments with target genome by calculating the distance of kmer frequecies
The script will split the target genome into fragments and select most different fragments with target genome by calculating the distance of kmer frequecies. Besides, the fragments overlapping with simple repeat sequences are removed using trf.
```
Usage: GIFEHGT kmerFilter [options] --fasta <genome_file>

Necessary input description:

  genome_file	<string>	The fasta file of whole target genome.

Options (defaults in parentheses):

    --help			Print this usage page.

  Output:

    --outdir	<string>	The result files will be output to this directory. (./kmerFilter)

  Other paramerters:

    --mode	<string>	standard or accelerated mode. (standard)

    --length	<int>		The length of split fragments. (1000)
     
    --step	<int>           The length of step to split the target genome. (800)

    --kmer	<int>		The length of kmer. (4)
     
    --kmerPer	<int>		The percentage of fragments selected that are most different from the target genome. (20)

    --trfPer	<int>           The percentage of sequences overlapping with simple repeats, and sequences above this percentage are judged as repeats and deleted as a whole. (50)
     
    --trfPara	<string>	The parameters of program trf (a comma-separated string), including matching weight,mismatching penaltyindel penalty,match probability,indel probability,minimum alignment score to report,maximum period size to report,maximum TR length expected. (2,5,7,80,10,50,1000,3)
```
Example command:
```
$ GIFEHGT kmerFilter --fasta path/to/example/data/fa/sacCer3_NC_001144.5.fa --mode accelerated
```
Results can be found in the directory `path/to/example/kmerFilter`.  
- `path/to/example/kmerFilter/4merFilter_removeSimRep.fa`: the fragments which are most different with target genome.
#### 2. `splitDB`: Split the genomic database into three groups: SG, CRG and DRG
The script will split the genomic database into three groups: SG, CRG and DRG and update taxonomic information of all organisms in genomic datasets.
```
Usage: GIFEHGT splitDB --id <genome_id> --taxo <genome_taxonomy> --taxoes <taxonomy_file> --type <type_file>

Necessary input description:

  genome_id		<string>        The genome id of target genome.

  genome_taxonomy	<string>	The taxonomy of target genome (a comma-separated string), whose format is like [kingdom,phylum,class,order].

  taxonomy_file         <string>        The taxonomy of all genomes. each line of which is like [genome_id	genome_kingdom 	genome_phylum	genome_class	genome_order].

  type_file		<string>	The kingdom of all genomes. Like [GCF_000001405.40_GRCh38.p14_genomic      Metazoa] for human.

Options (defaults in parentheses):

    --help                              Print this usage page.

  Output:

    --outdir		<string>	The result files will be output to this directory. (./splitDB)

  Other parameters:

    --distant		<string>	The distantly related group. kingdom or phylum can be chosen. (kingdom)
```
Example command:
```
$ GIFEHGT splitDB --id GCF_000146045.2 --taxo Fungi,Ascomycota,Saccharomycetes,Saccharomycetales --taxoes path/to/example/data/db/all_range.txt --type path/to/example/data/db/all_type.txt
```
Results can be found in the directory `path/to/example/splitDB`.  
#### 3. `seqAlign`: Sequence aligment using LASTZ
The script will acquire sequences of target genome which can align with genomes in DRG.
```
Usage: GIFEHGT seqAlign [options] --db <db_genome_dir>

Necessary input description:

  db_genome_dir		<string>	The directory of all genomes saved.

Options (defaults in parentheses):

    --help				Print this usage page.
  
  Input:

    --fragment		<string>	The fasta file of filtered fragments of target genome. (./kmerFilter/split.fasta)

    --genome		<string>	The fasta file of target genome. (./kmerFilter/new.fasta)

    --dbinfodir		<string>	The directory that all genome information saved. (./splitDB/)

  Output:

    --outdir		<string>	The result files will be output to this directory.

  Other parameters:

    --suffix		<string>	The suffix of genome file in database. (fna)

    --distant		<string>        The distantly related group. kingdom or phylum can be chosen. (kingdom)

    --identity		<float>		The threshold of identity. Two sequences above this threshold are considered to be similar. (0.5)
```
Example command:
```
$ GIFEHGT seqAlign --db path/to/example/data/db/refseq
```
Results can be found in the directory `path/to/example/seqAlign`.  
#### 4.  `screenHGT`: Screen potential HGTs using sequence alignment results
The script will acquire sequences of target genome which have higher identity with genomes in DRG than that in CRG. If Strict mode is chosen, it will require the similarity between sequences of target genome and genomes in CRG or DRG is high than one threshold. Besides, the sequences with too high or too low GC percentage, overlapped simple repeat, low complex repeat, single-copy-gene common to eukaryotes, ERV and mitochondrial and chloroplast (if have) sequences are removed. 
```
Usage: GIFEHGT screenHGT [options] --repeat <repeat_file> --singlecopy <singlecopy_file> --mitChl <mitChl_file> --taxo <genome_taxonomy>

Necessary input description:

  repeat_file		<string>	The repeat annotation file of target genome.

  singlecopy_file	<string>	The protein sequences file commom  to eukaryotes provided by BUSCO.

  mitChl_file		<string>	The fasta file of mitochondrial and chloroplast (if have) sequences.

  genome_taxonomy       <string>        The taxonomy of target genome (a comma-separated string), whose format is like [kingdom,phylum,class,order].
  
Options (defaults in parentheses):

    --help				Print this usage page.
 
  Input:

    --genome		<string>        The fasta file of target genome. (./kmerFilter/new.fasta)

    --dbInfoDir		<string> 	The directory that all genome information saved.

    --seqAlignDir	<string>        The directory that sequence alignment results saved. (./seqAlign/)

    --type		<string>        The type file updated by the taxanomy of target genome. (./splitDB/type.txt)

  Output:

    --outdir		<string>        The result files will be output to this directory. (./screenHGT/)

  Other parameters:

    --mode		<string>	Filter mode. Original or Strict can be chosen. (Original)

    --distant		<string>        The distantly related group. kingdom or phylum can be chosen. (kingdom)

    --self		<string>	The self group. phylum, classs, order or species can be chosen. (all of them)

    --length		<int>		The minimum length of HGTs. (135)

    --coverage		<float>		The minimum coverage of similar CRG sequences and DRG sequences. (0.6)

    --idenCdHitEst	<float>		The identity threshold for cd-hit-est. (0.8)
  
    --simCRG		<float>		The lowest similarity between homologous sequences in CRG and HGTs in Strict mode. (0.5)

    --simDRG		<float>		The lowest similarity between homologous sequences in DRG and HGTs in Strict mode. (0.6)
```
Example command (strict mode is used):
```
$ GIFEHGT screenHGT --mode Strict --repeat path/to/example/data/rmsk/rmsk.txt --singlecopy path/to/example/data/singlecopy/eukaryota.faa --mitChl path/to/example/data/mitochondria_chloroplast/mitochondria.fa --taxo Fungi,Ascomycota,Saccharomycetes,Saccharomycetales
```
Results can be found in the directory `path/to/example/screenHGT`. 
- `path/to/example/screenHGT/modeStrict/HGT.fa`: the sequences of potential HGTs.
#### 5.  `WGSValidate`: Validate potential HGTs using WGS datasets
The script will validate potential HGTs using WGS datasets.
```
Usage: GIFEHGT WGSValidate [options] --dbdir <db_genome_dir> --dbInfo <db_info_file> --taxo <genome_taxonomy> --dbWGSDir <db_WGSdata_dir>

Necessary input description:

  db_genome_dir			<string>	The directory of all genome database.

  db_info_file          	<string>        The information file of genomes in database.

  genome_taxonomy       	<string>        The taxonomy of target genome (a comma-separated string), whose format is like [kingdom,phylum,class,order].

  db_WGSdata_dir		<string>	The directory of WGS datasets.

Options (defaults in parentheses):

    --help                                      Print this usage page.

  Input:

    --genome			<string>        The fasta file of target genome. (./kmerFilter/new.fasta)

    --HGTId			<string>        The file of potential HGTs id. (./screenHGT/HGT.id or ./screenHGT/modeStrict/HGT.id)

    --HGTInfoDir		<string>        The information directory of homologous sequences for HGTs. (./screenHGT/tree_kingdom)

    --HGTInfoFile		<string>        The information file of homologous sequences for HGTs. (./screenHGT/HGT.info or ./screenHGT/modeStrict/HGT.info)

    --type			<string>        The type file of genomes in database updated before. (./splitDB/type.txt)

    --dbIdDir			<string>        The directory of genome id in database updated before. (./splitDB/)

  Output:

    --outdir			<string>        The result files will be output to this directory. (./WGSValidation/)

  Other paramerters:

    --mode                  	<string>        Filter mode. Original or Strict can be chosen. (Original)

    --length			<int>		The length HGT upstream and downstream sequence verified by WGS datasets. (150)

    --depth			<int>		The lowest depth at which the HGTs is covered by WGS datasets. (10)

    --distant			<string>	The distantly related group. kingdom or phylum can be chosen. (kingdom)

    --self			<string>        The self group. phylum, classs, order or species can be chosen. (all of them)

    --idenCdHitEst		<float>         The identity threshold for cd-hit-est. (0.8)
```
Example command (strict mode is used):
```
$ GIFEHGT WGSValidate --dbdir path/to/example/data/db/refseq/ --mode Strict --dbInfo path/to/example/data/db/all_info.txt --taxo Fungi,Ascomycota,Saccharomycetes,Saccharomycetales --dbWGSDir path/to/example/data/db/WGSdata/
```
Results can be found in the directory `path/to/example/WGSValidate`. 
- `path/to/example/WGSValidate/afterWGS/HGT.fa`: the sequences of final HGTs.
#### 6. `conPhyTree`: Construct sequence phylogenetic tree to validate HGTs
The script will construct sequence phylogenetic tree to validate HGTs.
```
Usage: GIFEHGT conPhyTree --genomeId <genome_id> --fullName <target_fullname> --dbId <db_id_file> --dbInfo <db_info_file>

Necessary input description:

  genome_id			<string>	The genome id of target genome.

  target_fullname		<string>	The full name of target species.

  db_id_file			<string>	The id of all genomes in database.

  db_info_file			<string>	The id and name of all genomes in database.

Options (defaults in parentheses):

    --help                                  	Print this usage page.

  Input:

    --HGTId			<string>	The file of potential HGTs id. (./WGSValidation/afterWGS/HGT.id)

    --HGTFa			<string>	The fasta file of potential HGTs. (./WGSValidation/afterWGS/HGT.fa)

    --HGTInfoDir		<string>        The information directory of homologous sequences for HGTs. (./screenHGT/tree/)

    --hitFile			<string>	The file name of homologous sequences for HGTs. (cov-hit.species.WGS.txt or cov-hit.species.strict.WGS.txt)

  Output:

    --outdir			<string>	The result files will be output to this directory. (./conPhyTree/)

  Other parameters:
    
    --mode              	<string>        Filter mode. Original or Strict can be chosen. (Original)

    --distant                  	<string>	The distantly related group. kingdom or phylum can be chosen. (kingdom)
```
Example command:
```
$ GIFEHGT conPhyTree --mode Strict --genomeId GCF_000146045.2 --fullName Saccharomyces_cerevisiae_S288C --dbId path/to/example/data/db/all_id.txt --dbInfo path/to/example/data/db/all_info.txt
```
Results can be found in the directory `path/to/example/conPhyTree`. 

[1]: http://www.bx.psu.edu/~rsharris/lastz
[2]: https://tandem.bu.edu/trf/downloads
[3]: https://github.com/weizhongli/cdhit/releases
[4]: https://bowtie-bio.sourceforge.net/bowtie2/index.shtml
[5]: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
[6]: https://github.com/samtools/samtools/releases/
[7]: https://mafft.cbrc.jp/alignment/software/
[8]: https://github.com/inab/trimal/releases
[9]: http://www.iqtree.org/
