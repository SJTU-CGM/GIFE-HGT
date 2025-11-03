#!/usr/bin/perl
package conPhyTree;

sub conPhyTree{
use strict;
use warnings;
use Getopt::Long;

my $usage="\nUsage: GIFEHGT conPhyTree --genomeId <genome_id> --fullName <target_fullname> --dbId <db_id_file> --dbInfo <db_info_file>

conPhyTree program is used to construct sequence phylogenetic tree to validate HGTs.

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

";

#Read screenHGT parameters
my $HGT_id_file = "./WGSValidation/afterWGS/HGT.id";
my $HGT_fa_file = "./WGSValidation/afterWGS/HGT.fa";
my $HGT_homologous_info_dir = "./screenHGT/tree/";
my $hit_file = "cov-hit.species.WGS.txt";
my $out_dir = "./conPhyTree/";
my $mode = "Original";
my $distant = "kingdom";
my $genome_id;
my $target_fullname;
my $db_id_file;
my $db_info_file;
my $help;
GetOptions(
    'HGTId=s'		=> \$HGT_id_file,
    'HGTFa=s'		=> \$HGT_fa_file,
    'HGTInfoDir=s'	=> \$HGT_homologous_info_dir,
    'hitFile=s'		=> \$hit_file,
    'genomeId=s'	=> \$genome_id,
    'fullName=s'	=> \$target_fullname,
    'dbId=s'		=> \$db_id_file,
    'dbInfo=s'		=> \$db_info_file,
    'help'      	=> \$help,
    'outdir=s'		=> \$out_dir,
    'mode=s'            => \$mode,
    'distant=s'		=> \$distant
) or die $usage;
die $usage if defined($help);
die $usage if(!$genome_id);
die $usage if(!$target_fullname);
die $usage if(!$db_id_file);
die $usage if(!$db_info_file);
die $usage if($mode ne "Original" && $mode ne "Strict");

if($mode eq "Strict"){
    $hit_file = "cov-hit.species.strict.WGS.txt";
}

#Check existence of output directory
if(-e $out_dir){
    print("Warning: output directory \"$out_dir\" already exists. Existing files will be overwritten.\n");
    `rm -rf $out_dir`;
}

#Adjust directory names and create output directory
unless($out_dir=~/\/$/){
    $out_dir.="/";
}
unless($HGT_homologous_info_dir=~/\/$/){
    $HGT_homologous_info_dir.="/";
}
my $out_dir_genomeid = $out_dir."tree_genomeid/";
my $out_dir_species = $out_dir."tree_species/";
mkdir($out_dir);
mkdir($out_dir_genomeid);
mkdir($out_dir_species);

#Acquire Homologous sequences of HGTs in all species
my %fasta = ();
my $fatitle = "";
open(FA,$HGT_fa_file)||die("conPhyTree: error with opening $HGT_fa_file\n");
while(<FA>){
    chomp();
    if($_ =~ />(.*)/){
        $fatitle = $1;
    }else{
	if(exists($fasta{$fatitle})){
	    $fasta{$fatitle} = $fasta{$fatitle}.$_;
	}else{
	    $fasta{$fatitle} = $_;
	}
    }
}
close FA;
open(ID,$HGT_id_file)||die("conPhyTree: error with opening $HGT_id_file\n");
while(<ID>){
    chomp();
    my $id = $_;
    my $out_dir_id = $out_dir_genomeid.$id."/";
    mkdir($out_dir_id);
    my $out_model_fa = $out_dir_id.$target_fullname.".fa";
    my $out_all_fa = $out_dir_id."all.fa";
    open(MFA,">$out_model_fa")||die("conPhyTree: error with writing to $out_model_fa\n");
    open(ALLFA,">$out_all_fa")||die("conPhyTree: error with writing to $out_all_fa\n");
    print MFA ">$target_fullname\n$fasta{$id}\n";
    print ALLFA ">$target_fullname\n$fasta{$id}\n";
    my $hit_file_id = $HGT_homologous_info_dir.$id."/".$hit_file;
    open(HIT,$hit_file_id)||die("conPhyTree: error with opening $hit_file_id\n");
    while(<HIT>){
	chomp();
	my @arr = split(/\t/,$_);
	if($arr[0] !~ $genome_id && $arr[0] ne "hit_species"){
	    print ALLFA ">$arr[0]\n";
	    my $hitfa = $HGT_homologous_info_dir.$id."/hit/".$arr[0].".fa";
            open(HITFA,$hitfa)||die("conPhyTree: error with opening $hitfa\n");
	    while(<HITFA>){
	        if($_ !~ />/){
		    print ALLFA "$_";
	        }
	    }
	}
    }
    close MFA; close ALLFA; close HIT; close HITFA;
    my $mafft = $out_dir_id."all.mafft";
    my $trimal = $out_dir_id."all.mafft.trimal";
    #multiple alignment
    my $commafft = "mafft --thread 4 --auto $out_all_fa > $mafft";
    #trimming
    my $comtrimal = "trimal -in $mafft -out $trimal -automated1";
    #construct phylogenetic tree
    my $comiqtree = "iqtree --threads-max 4 -st DNA -s $trimal -alrt 1000 -bb 1000";
    system($commafft);
    system($comtrimal);
    system($comiqtree);
}

#Covert genome id to species name
convertId2Name($db_id_file,$db_info_file,$out_dir_genomeid,$out_dir_species);
}

#Covert genome id to species name
sub convertId2Name{
#my ($target_abb,$target_fullname,$db_id_file,$db_info_file,$out_dir_genomeid,$out_dir_species) = @_;
my ($db_id_file,$db_info_file,$out_dir_genomeid,$out_dir_species) = @_;
open(DBID,$db_id_file)||die("conPhyTree.convertId2Name: error with opening $db_id_file\n");
my %id = ();
while(<DBID>){
    chomp();
    my @arr = split(/_/,$_);
    $id{"$arr[0]_$arr[1]"} = $_;
}
open(ID2SPECIES,$db_info_file)||die("conPhyTree.convertId2Name: error with opening $db_info_file\n");
my %hash = ();
while(<ID2SPECIES>){
    chomp();
    my @arr = split(/;/,$_);
    $hash{$id{$arr[0]}} = $arr[1];
}
opendir(DIR,$out_dir_genomeid) or die "conPhyTree.convertId2Name: error with opening $out_dir_genomeid\n";
my @dir = readdir DIR;
for my $d (@dir){
    if(!($d=~/.txt/ || $d=~/.tree/ || $d eq "." || $d eq "..")){
        opendir(NDIR,"$out_dir_genomeid$d") or die "conPhyTree.convertId2Name: error with opening $out_dir_genomeid$d";
        my @ndir = readdir NDIR;
        for my $file (@ndir){
            if($file=~/.treefile/){
                $out=$out_dir_species.$d.".tree";
                open(OUT,">$out")||die("conPhyTree.convertId2Name: error with writing to $out\n");
                open(ID,"$out_dir_genomeid$d/$file")||die("conPhyTree.convertId2Name: error with opening $out_dir_genomeid$d/$file\n");
                while(<ID>){
                    chomp();
                    $_=~s/(GC._.*?):/$hash{$1}:/g;
                    print OUT "$_\n";
                }
            }
        }
    }
}
close DBID; close ID2SPECIES; close DIR; close NDIR; close OUT;
}
1;
