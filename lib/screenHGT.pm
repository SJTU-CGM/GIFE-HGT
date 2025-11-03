#!/usr/bin/perl
package screenHGT;

sub screenHGT{
use strict;
use warnings;
use lib '.';
use base;
use Getopt::Long;

my $usage="\nUsage: GIFEHGT screenHGT [options] --repeat <repeat_file> --singlecopy <singlecopy_file> --mitChl <mitChl_file> --taxo <genome_taxonomy>

screenHGT program is used to screen potential HGTs from sequence alignment results.

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

    --self		<string>	The self group. phylum, class, order or species can be chosen. (all of them)

    --length		<int>		The minimum length of HGTs. (135)

    --coverage		<float>		The minimum coverage of similar CRG sequences and DRG sequences. (0.6)

    --idenCdHitEst	<float>		The identity threshold for cd-hit-est. (0.8)
  
    --simCRG		<float>		The lowest similarity between homologous sequences in CRG and HGTs in Strict mode. (0.5)

    --simDRG		<float>		The lowest similarity between homologous sequences in DRG and HGTs in Strict mode. (0.6)

";

#Read screenHGT parameters
my $genome_file = "./kmerFilter/new.fasta";
my $seq_alignment_dir = "./seqAlign/";
my $type_file = "./splitDB/type.txt";
my $out_dir = "./screenHGT/";
my $screen_distant = "kingdom";
my $screen_self = "phylum,class,order,species";
my $screen_length = 135;
my $screen_coverage = 0.6;
my $cdhit_threshold = 0.8;
my $sim_close = 0.5;
my $sim_distant = 0.6;
my $mode = "Original";
my $db_dir = "./splitDB/";
my $repeat_file;
my $singlecopy_file;
my $mitChl_file;
my $taxanomy_target;
my $help;
GetOptions(
    'genome=s'		=> \$genome_file,
    'seqAlign=s'	=> \$seq_alignment_dir,
    'dbInfoDir=s'	=> \$db_dir,
    'repeat=s' 		=> \$repeat_file,
    'singlecopy=s'	=> \$singlecopy_file,
    'mitChl=s'		=> \$mitChl_file,
    'type=s' 		=> \$type_file,
    'taxo=s'		=> \$taxanomy_target,
    'outdir=s'  	=> \$out_dir,
    'help'      	=> \$help,
    'mode=s'      	=> \$mode,
    'distant=s' 	=> \$screen_distant,
    'self=s'		=> \$screen_self,
    'length=i'		=> \$screen_length,
    'coverage=f'	=> \$screen_coverage,
    'idenCdHitEst=f'	=> \$cdhit_threshold,
    'simCRG=f'		=> \$sim_close,
    'simDRG=f'  	=> \$sim_distant
) or die $usage;

die $usage if defined($help);
die $usage if(!$db_dir);
die $usage if(!$repeat_file);
die $usage if(!$singlecopy_file);
die $usage if(!$mitChl_file);
die $usage if(!$taxanomy_target);
die $usage if($mode ne "Original" && $mode ne "Strict");

#Check existence of output directory
if(-e $out_dir){
    print("Warning: output directory \"$out_dir\" already exists. Existing files will be overwritten.\n");
    `rm -rf $out_dir`;
}

#Adjust directory names and create directory
unless($seq_alignment_dir=~/\/$/){
    $seq_alignment_dir.="/";
}
unless($db_dir=~/\/$/){
    $db_dir.="/";
}
unless($out_dir=~/\/$/){
    $out_dir.="/";
}
mkdir($out_dir);
my $out_dir_close = $out_dir."close/";
my $out_dir_distant = $out_dir."distant/";
my $out_dir_singlecopy = $out_dir."singlecopy/";
my $out_dir_ERV = $out_dir."ERV/";
my $out_dir_mitChl = $out_dir."mitChl/";
my $out_dir_tree = $out_dir."tree/";
mkdir($out_dir_close);
mkdir($out_dir_distant);
mkdir($out_dir_singlecopy);
mkdir($out_dir_ERV);
mkdir($out_dir_mitChl);
mkdir($out_dir_tree);
my $closeid = $db_dir.$screen_distant.".id";
my $distantid = $db_dir."e_".$screen_distant.".id";
my ($genome_kingdom,$genome_phylum,$genome_class,$genome_order) = split(/,/,$taxanomy_target);

#Intersection of CRG and DRG similar sequences
intersect($genome_file,$seq_alignment_dir,$closeid,$distantid,$out_dir,$out_dir_close,$out_dir_distant,$screen_length,$screen_coverage);

#Remove NNN sequences
removeNNN($out_dir);

#Remove sequence with too high or too low GC percentage 
#removeGC($out_dir);

#Remove sequence half overlapped Simple Repeat and Low complex Repeat
removeRepeat($genome_file,$repeat_file,$out_dir);

#Remove single-copy-gene common to eukaryotes
removeSCG($singlecopy_file,$out_dir);

#Remove ERV
removeERV($genome_file,$repeat_file,$out_dir);

#Remove mitochondrial and chloroplast (if have) sequences
removeMitChl($mitChl_file,$out_dir);

#Get homology information of candidate HGT sequences in other species
my $noNnogcnormsknoSCGnoERVnoMC_bed = $out_dir."noNnormsknoSCGnoERVnoMC.bed";
open(ALLID,$noNnogcnormsknoSCGnoERVnoMC_bed)||die("screenHGT: error with opening $noNnogcnormsknoSCGnoERVnoMC_bed\n");
while(<ALLID>){
    chomp();
    my @arr=split(/\t/,$_);
    my $id = "$arr[0]-$arr[1]-$arr[2]";
    getHomoInfo($seq_alignment_dir,$closeid,$id,$out_dir_tree);
}
close ALLID;

#Calculate species distribution of candidate HGT sequences
my $hit_file = "cov-hit.species.txt";
calDistribution($out_dir,$screen_distant,$genome_order,$genome_class,$genome_phylum,$genome_kingdom,$noNnogcnormsknoSCGnoERVnoMC_bed,$type_file,$out_dir_tree,$hit_file);

#The minimum identity of HGT in DRG is higher than the maximum identity in CRG
distantIdenHclose($screen_distant,$screen_self,$hit_file,$db_dir,$type_file,$out_dir_tree,$genome_file,$out_dir);

#Remove redundancy using cd-hit-est
my $dHcfa = $out_dir."distantHclose_".$screen_distant.".fa";
my $streeinfo = $out_dir."tree_".$screen_distant.".select.info";
rmRedundancy($dHcfa,$cdhit_threshold,$streeinfo,$out_dir);

#Strict mode is chosen
if($mode eq "Strict"){
my $out_dir_strict = $out_dir."modeStrict/";
mkdir($out_dir_strict);
my $HGTbed = $out_dir."HGT.bed";
open(HGTBED,$HGTbed)||die("screenHGT.modeStrict: error with opening $HGTbed\n");
while(<HGTBED>){
    chomp();
    my @arr=split(/\t/,$_);
    my $id = "$arr[0]-$arr[1]-$arr[2]";
    screenHomoInfo($closeid,$id,$out_dir_tree,$sim_close,$sim_distant);
}
$hit_file = "cov-hit.species.strict.txt";
calDistribution($out_dir_strict,$screen_distant,$genome_order,$genome_class,$genome_phylum,$genome_kingdom,$HGTbed,$type_file,$out_dir_tree,$hit_file);
distantIdenHclose($screen_distant,$screen_self,$hit_file,$db_dir,$type_file,$out_dir_tree,$genome_file,$out_dir_strict);
my $dHcfa = $out_dir_strict."distantHclose_".$screen_distant.".fa";
my $streeinfo = $out_dir_strict."tree_".$screen_distant.".select.info";
rmRedundancy($dHcfa,$cdhit_threshold,$streeinfo,$out_dir_strict);
}

}

#Intersection of CRG and DRG similar sequences
sub intersect{
use lib '.';
use base;
my ($genome_file,$seq_alignment_dir,$closeid,$distantid,$out_dir,$out_dir_close,$out_dir_distant,$screen_length,$screen_coverage) = @_;

open(CID,$closeid)||die("screenHGT.intersect: error with opening $closeid\n");
open(DID,$distantid)||die("screenHGT.intersect: error with opening $distantid\n");

#merge CRG
while(<CID>){
    chomp();
    my $file = $seq_alignment_dir."all/iden/".$_.".bed";
    my $file_sort = $seq_alignment_dir."all/iden/".$_.".bed.sort";
    my $file_merge = $out_dir_close.$_.".bed";
    base::sort($file,$file_sort);
    base::merge($file_sort,$file_merge);
}
close CID;

#merge DRG
while(<DID>){
    chomp();
    my $file = $seq_alignment_dir."all/iden/".$_.".bed";
    my $file_sort = $seq_alignment_dir."all/iden/".$_.".bed.sort";
    my $file_tmp = $out_dir_distant."tmp.txt";
    my $file_merge = $out_dir_distant.$_.".bed";
    base::sort($file,$file_sort);
    base::merge($file_sort,$file_tmp);
    open(TMP,$file_tmp)||die("screenHGT.intersect: error with opening tmp.txt\n");
    open(MERGE,">$file_merge")||die("screenHGT.intersect: error with opening $file_merge\n");
    while(<TMP>){
	chomp();
	my @arr = split(/\t/,$_);
	if($arr[2]-$arr[1]+1 >= $screen_length){
	    print MERGE "$_\n";
        }
    }
}
close DID; close TMP; close MERGE;

my $file_mergeall = $out_dir_distant."distant.txt";
my $comcat = "cat $out_dir_distant/*bed |sort -k1,1 -k2n,2 | uniq > $file_mergeall";
system($comcat);

#intersect CRG and DRG
my $file_all = $out_dir."all.out";
open(DISTANT,$file_mergeall)||die("screenHGT.intersect: error with opening $file_mergeall\n");
open(OUT,">$file_all")||die("screenHGT.intersect: error with writing $file_all\n");
opendir(CDIR,$out_dir_close)||die("screenHGT.intersect: error with opening $out_dir_close\n");

while(<DISTANT>){
    chomp();
    my ($chr_this,$start_this,$end_this) = split(/\s+/,$_);
    my $length_this = $end_this-$start_this+1;
    my $cov_species = 0;

    while(<CDIR>){
        chomp();
        if($_ =~ /.*bed/){
      	    my $file = $_;
        }
        open(FILE,$file)||die("error with opening $file\n");
        my $cov_length = 0;

        while(<FILE>){
            chomp();
            my ($chr,$start,$end) = split(/\s+/,$_);
            if($chr gt $chr_this){
                last;
            }elsif($chr eq $chr_this){
                if($start >= $end_this){
                    last;
                }elsif(!($end <= $start_this)){
                    my ($cov_start,$cov_end) = ($start_this,$end_this);
                    if($start > $cov_start){$cov_start = $start;}
                    if($end < $cov_end){$cov_end = $end;}
                    $cov_length += $cov_end-$cov_start+1;
                }
            }
        }
        if($cov_length/$length_this >= $screen_coverage){
            $cov_species += 1;
        }
        close FILE;
    }
    print OUT "$chr_this\t$start_this\t$end_this\t$cov_species\n";
}

close DISTANT; close OUT; close CDIR;

#Merge bed
my $file_allbed = $out_dir."all.bed.mid";
my $file_allmerge = $out_dir."all.bed";
base::id2bed($file_all,$file_allbed);
base::merge($file_allbed,$file_allmerge);
my $file_allfasta = $out_dir."all.fa";
base::getSeq($genome_file,$file_allmerge,$file_allfasta);
}

#Remove sequences with NNN
sub removeNNN{
use lib '.';
use base;
my ($out_dir) = @_;
my $allfasta = $out_dir."all.fa";
my $noNfasta = $out_dir."noN.fa";
my $noNbed = $out_dir."noN.bed";

open(ALL,$allfasta)||die("screenHGT.removeNNN: error with opening $allfasta\n");
open(NON,">$noNfasta")||die("screenHGT.removeNNN: error with writing $noNfasta\n");
my $id = "";
while(<ALL>){
    chomp();
    if($_=~ />/){
        $id = $_;
    }
    elsif(!($_ =~ /N/)){
        print NON "$id\n$_\n";
    }
}
close ALL; close NON;

base::fastatobed($noNfasta,$noNbed);
}

#Remove sequence with too high or too low GC percentage
#sub removeGC{
#use lib '.';
#use base;
#use kmerFilter;
#my ($out_dir) = @_;
#my $noNfasta = $out_dir."noN.fa";
#my $noNfa_gc = $out_dir."noN-gc.txt";
#my $noNnogc_bed = $out_dir."noNnogc.bed";

#my @kmer = ('A','T','C','G');
#kmerFilter::kmerFragments($noNfasta,1,$noNfa_gc,\@kmer);

#open(ALL,$noNfa_gc)||die("screenHGT.removeGC: error with opening $allfasta\n");
#open(OUT,">$noNnogc_bed")||die("screenHGT.removeGC: error with writing $noNnogc_bed\n");
#while(<ALL>){
#    chomp();
#    @arr = split(/\t|-|:/,$_);
#    if($arr[5]+$arr[6] < 0.7 && $arr[5]+$arr[6] > 0.3){
#        print OUT "$arr[0]\t$arr[1]\t$arr[2]\n";
#    }
#}
#close ALL; close OUT;
#}

#Remove sequence half overlapped Simple Repeat and Low complex Repeat
sub removeRepeat{
use lib '.';
use base;
#use removeSimRep;
my ($genome_file,$repeat_file,$out_dir) = @_;
my $noNnogc_bed = $out_dir."noN.bed";
my $noNnogc_rmsk = $out_dir."noN_rmsk.txt";
my $noNnogcnormskfilter_bed = $out_dir."noNnormskmid.filter.bed";
my $noNnogcnormskmid_bed = $out_dir."noNnormskmid.bed";
my $noNnogcnormskmid_fa = $out_dir."noNnormskmid.fa";
my $noNnogcnormskmid_faname = "noNnormskmid.fa";
my $noNnogcnormsk_name = "noNnormsk";
my $noNnogcnormsk_bed = $out_dir."noNnormsk.bed";
my $noNnogcnormsk_fa = $out_dir."noNnormsk.fa";

my $repeat_bed = $repeat_file.".bed";
open(REPEAT,$repeat_file)||die("screenHGT.removeRepeat: error with opening $repeat_file\n");
open(REPBED,">$repeat_bed")||die("screenHGT.removeRepeat: error with writing $repeat_bed\n");
while(<REPEAT>){
    chomp();
    my @arr = split(/\s+/,$_);
    if($arr[5] ne "begin" && $arr[0] ne ""){
    if($arr[8] eq "C" || $arr[8] eq "+"){
        print REPBED "$arr[4]\t$arr[5]\t$arr[6]\t$arr[10]\n";
    }else{
        print REPBED "$arr[5]\t$arr[6]\t$arr[7]\t$arr[11]\n";
    }
    }
}
close REPBED; close REPEAT;
base::biodiff($repeat_bed,$noNnogc_bed,$noNnogc_rmsk);
base::filterSimple($noNnogc_rmsk,$noNnogcnormskfilter_bed);
base::notin($noNnogcnormskfilter_bed,$noNnogc_bed,$noNnogcnormskmid_bed);
base::getSeq($genome_file,$noNnogcnormskmid_bed,$noNnogcnormskmid_fa);

#removeSimRep::removeSimRep($noNnogcnormskmid_fa,$noNnogcnormskmid_faname,$genome_file,$out_dir,$noNnogcnormsk_name);
rmSrpTrf($noNnogcnormskmid_fa,$noNnogcnormskmid_faname,$genome_file,$out_dir,$noNnogcnormsk_name);
my $tmpbed = $out_dir.$noNnogcnormsk_name."_removeSimRep.bed";
my $tmpfa = $out_dir.$noNnogcnormsk_name."_removeSimRep.fa";
my $commvbed = "mv $tmpbed $noNnogcnormsk_bed";
my $commvfa = "mv $tmpfa $noNnogcnormsk_fa";
system($commvbed);
system($commvfa);
}

#Remove Simple Repeat using trf strandard parameters
sub rmSrpTrf{
use lib '.';
use kmerFilter;
my ($kmer_filter_file,$kmer_filter_name,$genome_file,$out_dir,$trf_out_name) = @_;
my $trf_out_file = $out_dir.$trf_out_name.".trf.mask";
my $outbed = $out_dir.$trf_out_name."_removeSimRep.bed";
my $outfa = $out_dir.$trf_out_name."_removeSimRep.fa";
my $remove_per = 50;
my $remove_match = 2;
my $remove_Mismatch = 5;
my $remove_Delta = 7;
my $remove_PM = 80;
my $remove_PI = 10;
my $remove_Minscore = 50;
my $remove_MaxPeriod = 1000;
my $remove_len = 3;
#Run trf to mask the simple repeat sequences
kmerFilter::runtrf($kmer_filter_file,$kmer_filter_name,$trf_out_file,$remove_match,$remove_Mismatch,$remove_Delta,$remove_PM,$remove_PI,$remove_Minscore,$remove_MaxPeriod,$remove_len);
#Remove fragments with simple repeat sequences
kmerFilter::removeSim($trf_out_file,$remove_per,$outbed,$genome_file,$outfa);
}

#Remove single-copy-gene common to eukaryotes
sub removeSCG{
use lib '.';
use base;
my ($singlecopy_file,$out_dir) = @_;
my $db_eukaryota = $out_dir."singlecopy/eukaryota";
my $noNnogcnormsk_fa = $out_dir."noNnormsk.fa";
my $noNnogcnormsk_bed = $out_dir."noNnormsk.bed";
my $out_blastx = $out_dir."singlecopy/HGT.blastx";
my $out_bed = $out_dir."singlecopy/HGT.singlecopy.bed";
my $noNnogcnormsknoSCG_bed = $out_dir."noNnormsknoSCG.bed";

my $commkdb = "makeblastdb -dbtype prot -in $singlecopy_file -out $db_eukaryota";
my $comblastx = "blastx -query $noNnogcnormsk_fa -db $db_eukaryota -out $out_blastx -evalue 1e-5 -outfmt 7 -num_threads 4";
system($commkdb);
system($comblastx);
open(ALL,$out_blastx)||die("screenHGT.removeSCG: error with opening $out_blastx\n");
open(OUT,">$out_bed")||die("screenHGT.removeSCG: error with writing $out_bed\n");
my %hash = ();
while(<ALL>){
    chomp();
    @arr = split(/\t/,$_);
    if($_ !~ /#/ && $arr[2]>=30 && not exists($hash{$arr[0]})){
	print "$arr[0]\t$arr[2]\n";
	$hash{$arr[0]} = 1;
	@bed = split(/-/,$arr[0]);
        print OUT "$bed[0]\t$bed[1]\t$bed[2]\n";
    }
}
close ALL; close OUT;

base::notin($out_bed,$noNnogcnormsk_bed,$noNnogcnormsknoSCG_bed);
}

#Remove ERV
sub removeERV{
use lib '.';
use base;
my ($genome_file,$repeat_file,$out_dir) = @_;
my $ERV_bed = $out_dir."ERV/ERV.bed";
my $ERV_fa = $out_dir."ERV/ERV.fa";
my $ERV_db = $out_dir."ERV/ERV";
my $ERV_blast = $out_dir."ERV/blast.txt";
my $noNnogcnormsknoSCG_bed = $out_dir."noNnormsknoSCG.bed";
my $noNnogcnormsknoSCG_fa = $out_dir."noNnormsknoSCG.fa";
my $noNnogcnormsknoSCG_ERV = $out_dir."noNnormsknoSCG_ERV";
my $noNnogcnormsknoSCG_ERV_bed = $out_dir."noNnormsknoSCG_ERV.bed";
my $noNnogcnormsknoSCGnoERV_bed = $out_dir."noNnormsknoSCGnoERV.bed";
my $noNnogcnormsknoSCGnoERV_fa = $out_dir."noNnormsknoSCGnoERV.fa";
open(REPEAT,$repeat_file)||die("screenHGT.removeERV: error with opening $repeat_file\n");
open(ERV,">$ERV_bed")||die("screenHGT.removeERV: error with writing $ERV_bed\n");
while(<REPEAT>){
    chomp();
    my @arr = split(/\t/,$_);
    if($arr[10] =~ /ERV/){
	print ERV "$arr[4]\t$arr[5]\t$arr[6]\n";
    }elsif($arr[11] =~ /ERV/){
	print ERV "$arr[5]\t$arr[6]\t$arr[7]\n";
    }
}
close REPEAT; close ERV;
base::getSeq($genome_file,$ERV_bed,$ERV_fa);
base::getSeq($genome_file,$noNnogcnormsknoSCG_bed,$noNnogcnormsknoSCG_fa);

base::biodiff($ERV_bed,$noNnogcnormsknoSCG_bed,$noNnogcnormsknoSCG_ERV);
open(ERVINFO,"$noNnogcnormsknoSCG_ERV")||die("screenHGT.removeERV: error with opening $noNnogcnormsknoSCG_ERV\n");
open(ERVBED,">$noNnogcnormsknoSCG_ERV_bed")||die("screenHGT.removeERV: error with writing $noNnogcnormsknoSCG_ERV_bed\n");
while(<ERVINFO>){
    chomp();
    my @arr = split(/\t| |;/,$_);
    print ERVBED "$arr[0]\t$arr[1]\t$arr[2]\n";
}
close ERVINFO;

my $commkdb = "makeblastdb -dbtype nucl -in $ERV_fa -out $ERV_db";
my $comblastn = "blastn -task blastn -query $noNnogcnormsknoSCG_fa -db $ERV_db -out $ERV_blast -evalue 1e-5 -outfmt 7 -num_threads 4";
system($commkdb);
system($comblastn);
open(ALL,$ERV_blast)||die("screenHGT.removeERV: error with opening $ERV_blast\n");
my %hash = ();
while(<ALL>){
    chomp();
    @arr = split(/\t/,$_);
    if($_ !~ /#/ && $arr[2]>=80 && $arr[3]>=100 && not exists($hash{$arr[0]})){
        $hash{$arr[0]} = 1;
        print ERVBED "$arr[0]\n";
    }
}
close ALL; close ERVBED;
base::notin($noNnogcnormsknoSCG_ERV_bed,$noNnogcnormsknoSCG_bed,$noNnogcnormsknoSCGnoERV_bed);
base::getSeq($genome_file,$noNnogcnormsknoSCGnoERV_bed,$noNnogcnormsknoSCGnoERV_fa);
}

#Remove mitochondrial and chloroplast (if have) sequences
sub removeMitChl{
use lib '.';
use base;
my ($mitChl_file,$out_dir) = @_;
my $db_mitChl = $out_dir."mitChl/mitChl";
my $noNnogcnormsknoSCGnoERV_fa = $out_dir."noNnormsknoSCGnoERV.fa";
my $noNnogcnormsknoSCGnoERV_bed = $out_dir."noNnormsknoSCGnoERV.bed";
my $out_blast = $out_dir."mitChl/HGT.blast";
my $out_bed = $out_dir."mitChl/HGT.mitChl.bed";
my $noNnogcnormsknoSCGnoERVnoMC_bed = $out_dir."noNnormsknoSCGnoERVnoMC.bed";

my $commkdb = "makeblastdb -dbtype nucl -in $mitChl_file -out $db_mitChl";
my $comblast = "blastn -query $noNnogcnormsknoSCGnoERV_fa -db $db_mitChl -out $out_blast -evalue 1e-5 -outfmt 7 -num_threads 4";
system($commkdb);
system($comblast);
open(ALL,$out_blast)||die("screenHGT.removeMitChl: error with opening $out_blast\n");
open(OUT,">$out_bed")||die("screenHGT.removeMitChl: error with writing $out_bed\n");
my %hash = ();
while(<ALL>){
    chomp();
    @arr = split(/\t/,$_);
    if($_ !~ /#/ && $arr[2]>=80 && not exists($hash{$arr[0]})){
        $hash{$arr[0]} = 1;
        @bed = split(/-/,$arr[0]);
        print OUT "$bed[0]\t$bed[1]\t$bed[2]\n";
    }
}
close ALL; close OUT;

base::notin($out_bed,$noNnogcnormsknoSCGnoERV_bed,$noNnogcnormsknoSCGnoERVnoMC_bed);
}

#Get homology information of candidate HGT sequences in other species
sub getHomoInfo{
my ($seqAlign_dir,$closeid,$hgtid,$out_dir_tree) = @_;
my ($hgtchr,$hgtstart,$hgtend) = split(/-/,$hgtid);
my $seqAlign_dir_cov = $seqAlign_dir."all/cov/";
my $out_dir_tree_id = $out_dir_tree.$hgtid."/";
my $out_dir_tree_id_hit = $out_dir_tree.$hgtid."/hit/";
my $out_cov = $out_dir_tree_id."cov.txt";
my $out_cov_best = $out_dir_tree_id."cov-hit.txt";
my $out_cov_best_species = $out_dir_tree_id."cov-hit.species.txt";
mkdir($out_dir_tree_id);
mkdir($out_dir_tree_id_hit);

my %close = ();
open(CLOSE,$closeid)||die("screenHGT.getHomoInfo: error with opening $closeid\n");
while(<CLOSE>){
    chomp();
    $close{$_} = 1;
}
close CLOSE;
open(COV,">$out_cov")||die("screenHGT.getHomoInfo: error with writing to $out_cov\n");
opendir(SEQA,$seqAlign_dir_cov)||die("screenHGT.getHomoInfo: error with opening $seqAlign_dir_cov\n");
my @files = readdir(SEQA);
foreach $f (@files){
    chomp();
    if($f =~ /(.*).bed.merge/){
        my $file = $seqAlign_dir_cov.$f;
	my $fileid = $1;
        open(FILE,$file)||die("screenHGT.getHomoInfo: error with opening $file\n");
        while(<FILE>){
	chomp();
	my @arr = split(/\t/,$_);
	my $cov = 1;
	if($_ =~ /$hgtchr/ && $arr[1] < $hgtend && $arr[2] > $hgtstart && $arr[7] > 0.5){
	    my ($start_new,$end_new) = ($hgtstart,$hgtend);
            if($arr[1] > $hgtstart){
                $start_new = $arr[1];
            }
            if($arr[2] < $hgtend){
                $end_new = $arr[2];
            }
            $cov = ($end_new-$start_new+1)/($hgtend-$hgtstart+1);
            my $iden = sprintf"%.4f",$arr[7]*$cov;
	    if(exists($close{$fileid})){
		print COV "$fileid\t$_\t$iden\n";
	    }elsif($arr[7] >= $iden_distant && $cov >= $cov_distant){
		print COV "$fileid\t$_\t$iden\n";
            }
        }
        }
        close FILE;
    }
}
close COV; close SEQA;

#best match for each species
open(COV,$out_cov)||die("error with opening $out_cov\n");
open(COVBEST,">$out_cov_best")||die("error with writing to $out_cov_best\n");
my %hash = ();
while(<COV>){
    chomp();
    my @arr = split(/\t/,$_);
    my ($species,$chr_this,$start_this,$end_this,$chr,$start,$end,$string,$identity,$similarity) = ($arr[0],$arr[1],$arr[2],$arr[3],$arr[4],$arr[5],$arr[6],$arr[7],$arr[8],$arr[9]);
    my $len = (abs($end_this-$start_this)+1) * $identity;
    if(not exists($hash{$species})){
        $hash{$species} = $chr_this.";".$start_this.";".$end_this.";".$chr.";".$start.";".$end.";".$identity.";".$similarity.";".$len.";".$string;
    }
    else{
        my $old = $hash{$species};
        my ($chr_this_old,$start_this_old,$end_this_old,$chr_old,$start_old,$end_old,$identity_old,$similarity_old,$len_old,$string_old) = split(/;/,$old);
        if($len > $len_old){
            $hash{$species} = $chr_this.";".$start_this.";".$end_this.";".$chr.";".$start.";".$end.";".$identity.";".$similarity.";".$len.";".$string;
        }
    }
}
foreach my $key(keys %hash){
    my ($chr_this,$start_this,$end_this,$chr,$start,$end,$identity,$similarity,$len,$string) = split(/;/,$hash{$key});
    print COVBEST "$key\t$chr\t$start\t$end\t$string\t$identity\t$similarity\t$chr_this\t$start_this\t$end_this\n";
}
close COV;close COVBEST;

#get species information
open(COVBEST,$out_cov_best)||die("error with opening to $out_cov_best\n");
open(COVBESTSPE,">$out_cov_best_species")||die("error with writing to $out_cov_best_species\n");
print COVBESTSPE "hit_species\tchr\tstart\tend\tstrand\tidentity\tsimilarity\n";
while(<COVBEST>){
    chomp();
    my ($id,$chr,$start,$end,$strand,$iden,$sim,$chr_match,$start_match,$end_match) = split(/\s+/,$_);
    if($start_match < $hgtstart){
        if($end_match <= $hgtend){
            my $ratio = ($hgtend-$hgtstart+1)/($end_match-$start_match+1);
            my $len_new = int($ratio*($end-$start+1));
            $start_match = $hgtstart;
            $start = $end-$len_new+1;
        }
        else{
            my $ratio = ($hgtend-$hgtstart+1)/($end_match-$start_match+1);
            my $len_change_left = int(($hgtstart-$start_match+1)/($end_match-$start_match+1)*($end-$start+1));
            my $len_change_right = int(($end_match-$hgtend+1)/($end_match-$start_match+1)*($end-$start+1));
            $start_match = $hgtstart;
            $end_match = $hgtend;
            $start = $start + $len_change_left;
            $end = $end - $len_change_right;
        }
    }
    else{
        if($end_match > $hgtend){
            my $ratio = ($hgtend-$start_match+1)/($end_match-$start_match+1);
            my $len_new = int($ratio*($end-$start+1));
            $end_match = $end_hgt;
            $end = $start+$len_new-1;
        }
    }
    print COVBESTSPE "$id\t$chr\t$start\t$end\t$strand\t$iden\t$sim\n";
}
close COVBEST; close COVBESTSPE;
}

#Filter homologous sequences of HGTs in Strict mode 
sub screenHomoInfo{
my ($closeid,$id,$out_dir_tree,$sim_close,$sim_distant) = @_;
my $file = $out_dir_tree.$id."/cov-hit.species.txt";
my $strictfile = $out_dir_tree.$id."/cov-hit.species.strict.txt";
open(FILE,$file)||die("screenHGT.screenHomoInfo: error with opening $file\n");
open(OUT,">$strictfile")||die("screenHGT.screenHomoInfo: error with writing to $strictfile\n");
my %close = ();
open(CLOSE,$closeid)||die("screenHGT.screenHomoInfo: error with opening $closeid\n");
while(<CLOSE>){
    chomp();
    $close{$_} = 1;
}
while(<FILE>){
    chomp();
    my @arr = split(/\t/,$_);
    ### for close and remote species, the identity setting is different
    if($arr[0] ne "hit_species"){
    if(($arr[6] >= $sim_distant && not exists($close{$arr[0]}))){
            print OUT "$_\n";
    }elsif( $arr[6] >= $sim_close && exists($close{$arr[0]})){
            print OUT "$_\n";
    }
    }else{
            print OUT "$_\n";
    }
}
close FILE;close OUT; close CLOSE;
}

#Calculate species distribution of candidate HGT sequences
sub calDistribution{
my ($out_dir,$screen_distant,$genome_order,$genome_class,$genome_phylum,$genome_kingdom,$noNnogcnormsknoSCGnoERVnoMC_bed,$type_file,$out_dir_tree,$hit_file) = @_;
my $treeinfo = $out_dir."tree_".$screen_distant.".info";
my $streeinfo = $out_dir."tree_".$screen_distant.".select.info";
my $streeid = $out_dir."tree_".$screen_distant.".select.id";
my @infohead = ();
push(@infohead,$genome_order);
push(@infohead,$genome_class);
push(@infohead,$genome_phylum);
push(@infohead,$genome_kingdom);
my @allkingdom = ("Bacteria","Viruses","Metazoa","Protozoa","Viridiplantae","Fungi");
foreach my $k (@allkingdom){
    if($k ne $genome_kingdom){
        push(@infohead,$k);
    }
}
open(TREEINFO,">$treeinfo")||die("screenHGT: error with writing to $treeinfo\n");
open(STREEINFO,">$streeinfo")||die("screenHGT: error with writing to $streeinfo\n");
open(STREEID,">$streeid")||die("screenHGT: error with writing to $streeid\n");
print TREEINFO "id\t";
print STREEINFO "id\t";
foreach my $h (@infohead){
    print TREEINFO "$h\t";
    print STREEINFO "$h\t";
}
print TREEINFO "all_species\n";
print STREEINFO "all_species\n";

open(ALLID,$noNnogcnormsknoSCGnoERVnoMC_bed)||die("screenHGT: error with opening $noNnogcnormsknoSCGnoERVnoMC_bed\n");
while(<ALLID>){
    chomp();
    my @arr=split(/\t|-/,$_);
    my $id = "$arr[0]-$arr[1]-$arr[2]";
    my $return = calDis($id,$type_file,$screen_distant,$out_dir_tree,$hit_file,\@infohead);
    my @return = @{$return};
    my $ifselected = pop(@return);
    my $spedis = join("\t",@return);
    print TREEINFO "$id\t$spedis\n";
    if($ifselected){
        print STREEINFO "$id\t$spedis\n";
        print STREEID "$id\n";
    }
}
close ALLID; close STREEID; close TREEINFO; close STREEINFO;
}

sub calDis{
my ($id,$type_file,$screen_distant,$out_dir_tree,$hit_file,$infohead) = @_;
my @infohead = @{$infohead};
my $cov_best_species = $out_dir_tree.$id."/".$hit_file;
open(COVBESTSPE,$cov_best_species)||die("screenHGT.calDis: error with opening $cov_best_species\n");
open(TYPE,$type_file)||die("screenHGT.calDis: error with opening $type_file\n");

my %hash = ();
while(<TYPE>){
    chomp();
    my @arr = split(/\s+/,$_);
    $hash{$arr[0]} = $arr[1];
}

my @num = (0,0,0,0,0,0,0,0,0);
while(<COVBESTSPE>){
    chomp();
    my $count = 0;
    my @arr = split(/\s+/,$_);
    if(exists($hash{$arr[0]})){
	foreach $i (@infohead){
	    if($hash{$arr[0]} eq $i){
                $num[$count]++;
	    }
	    $count++;
        }
    }
}
my $all = 0;
foreach $n (@num){
    $all = $all + $n;
}
#push($all,@num);
my $ifright = 0;
if($num[6]*$num[7]*$num[8]==0){
if($screen_distant eq "kingdom"){
    if($all-$num[0]-$num[1]-$num[2]-$num[3]>0){
        $ifright = 1;
    }
}elsif($screen_distant eq "phylum"){
    if($all-$num[0]-$num[1]-$num[2]-$num[3]-$num[4]>0){
        $ifright = 1;
    }
}elsif($screen_distant eq "class"){
    if($all-$num[0]-$num[1]-$num[2]-$num[3]-$num[4]-$num[5]>0){
        $ifright = 1;
    }
}
}
my @return = (@num,$all,$ifright);
return \@return;
close COVBESTSPE; close TYPE;
}

#The minimum identity of HGT in DRG is higher than the maximum identity in CRG
sub distantIdenHclose{
my ($screen_distant,$screen_self,$hit_file,$db_dir,$type_file,$out_dir_tree,$genome_file,$out_dir) = @_;

my @HGT_self = split(/,/,$screen_self);
my $streeid = $out_dir."tree_".$screen_distant.".select.id";
my @self = ();
if($screen_distant eq "kingdom"){
    @self = ("phylum","class","order","species");
}elsif($screen_distant eq "phylum"){
    @self = ("class","order","species");
}elsif($screen_distant eq "class"){
    @self = ("order","species");
}
my $sizeself = @self;
my @selfid = ();
my @closeids = ();
my $distantid = $db_dir."e_".$screen_distant.".id";
for(my $i=0;$i<$sizeself;$i++){
    $selfid[$i] = $db_dir.$self[$i].".id";
    $closeids[$i] = $db_dir.$screen_distant."-".$self[$i].".id";
}
for(my $i=0;$i<$sizeself;$i++){
    distantHclose($streeid,$self[$i],$selfid[$i],$closeids[$i],$distantid,$screen_distant,$type_file,$out_dir_tree,$hit_file,$out_dir);
}
for(my $i=0;$i<$sizeself-1;$i++){
    my $highlevel = $out_dir."distantHclose_".$screen_distant."-".$self[$i].".bed";
    my $lowlevel = $out_dir."distantHclose_".$screen_distant."-".$self[$i+1].".bed";
    my $tmp = $out_dir."tmp.txt";
    base::notin($lowlevel,$highlevel,$tmp);
    my $commv = "mv $tmp $highlevel";
    system($commv);
}
my $sizeHGTself = @HGT_self;
my $dHcbed = $out_dir."distantHclose_".$screen_distant.".bed";
for(my $i=0;$i<$sizeHGTself;$i++){
    `cat $out_dir"distantHclose_"$screen_distant"-"$HGT_self[$i]".bed" >> $dHcbed`;
    print $out_dir."distantHclose_".$screen_distant."-".$HGT_self[$i].".bed\n";
}
my $dHcfa = $out_dir."distantHclose_".$screen_distant.".fa";
base::getSeq($genome_file,$dHcbed,$dHcfa);
}

sub distantHclose{
my ($streeid,$self,$selfid,$closeids,$distantid,$screen_distant,$type_file,$out_dir_tree,$hit_file,$out_dir) = @_;
my $out_info = $out_dir."distantHclose_".$screen_distant."-".$self.".info";
my $out_bed = $out_dir."distantHclose_".$screen_distant."-".$self.".bed";

my %type = ();
open(TYPE,$type_file)||die("screenHGT.distantHclose: error with opening $type_file\n");
while(<TYPE>){
    chomp();
    my @arr = split(/\s+/,$_);
    $type{$arr[0]} = $arr[1];
}
my %distant = ();
open(DISTANT,$distantid)||die("screenHGT.distantHclose: error with opening $distantid\n");
while(<DISTANT>){
    chomp();
    $distant{$_} = 1;
}
my %close = ();
open(CLOSE,$closeids)||die("screenHGT.distantHclose: error with opening $closeids\n");
while(<CLOSE>){
    chomp();
    $close{$_} = 1;
}

open(ID,$streeid)||die("screenHGT.distantHclose: error with opening $streeid\n");
open(OUTINFO,">$out_info")||die("screenHGT.distantHclose: error with writing to $outinfo\n");
open(OUTBED,">$out_bed")||die("screenHGT.distantHclose: error with writing to $outbed\n");
while(<ID>){
    chomp();
    my $id = $_;
    my $file = $out_dir_tree.$id."/".$hit_file;
    open(FILE,$file)||die("screenHGT.distantHclose: error with opening $file\n");
    my $close_best = 0;
    my $remote_best = 0;
    my $remote_worst = 2;
    my $remote_species = "";
    while(<FILE>){
        my @arr = split(/\s+/,$_);
        if(exists($close{$arr[0]}) && $arr[5]>$close_best){
            $close_best = $arr[5];
        }
        if(exists($distant{$arr[0]}) && $arr[5]>$remote_best){
            $remote_best = $arr[5];
            $remote_species = $arr[0];
        }
        if(exists($distant{$arr[0]}) && $arr[5]>0 && $arr[5]<$remote_worst){
            $remote_worst = $arr[5];
        }

    }
    if($remote_worst<=1 &&  $remote_worst-$close_best>=0){
	my $bed = $id;
	$bed =~ s/-/\t/g;
	print OUTBED "$bed\n";
        print OUTINFO "$id\t$close_best\t$remote_best\t$remote_species\t$type{$remote_species}\n";
    }
    close FILE;
}

close ID; close TYPE; close CLOSE; close DISTANT;
close OUTBED; close OUTINFO;
}

#Remove redundancy using cd-hit-est
sub rmRedundancy{
my ($dHcfa,$cdhit_threshold,$streeinfo,$out_dir) = @_;
my $HGTfa = $out_dir."HGT.fa";
my $HGTbed = $out_dir."HGT.bed";
my $HGTid = $out_dir."HGT.id";
my $HGTinfo = $out_dir."HGT.info";
my $comcdhit = "cd-hit-est -i $dHcfa -c $cdhit_threshold -o $HGTfa";
system($comcdhit);
base::fastatobed($HGTfa,$HGTbed);
base::bed2id($HGTbed,$HGTid);
base::id2info($HGTid,$streeinfo,$HGTinfo);

my @self = ();
if($screen_distant eq "kingdom"){
    @self = ("phylum","class","order","species");
}elsif($screen_distant eq "phylum"){
    @self = ("class","order","species");
}elsif($screen_distant eq "class"){
    @self = ("order","species");
}
my $sizeself = @self;
for(my $i=0;$i<$sizeself;$i++){
    my $oribed = $out_dir."distantHclose_".$screen_distant."-".$self[$i].".bed";
    my $newbed = $out_dir."HGT_".$screen_distant."-".$self[$i].".bed";
    my $newid = $out_dir."HGT_".$screen_distant."-".$self[$i].".id";
    base::in($HGTbed,$oribed,$newbed);
    base::bed2id($newbed,$newid);
}
}

1;
