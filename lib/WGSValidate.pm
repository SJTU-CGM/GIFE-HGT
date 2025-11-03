#!/usr/bin/perl
package WGSValidate;

sub WGSValidate{
use strict;
use warnings;
use lib '.';
use screenHGT;
use Getopt::Long;

my $usage="\nUsage: GIFEHGT WGSValidate [options] --dbdir <db_genome_dir> --dbInfo <db_info_file> --taxo <genome_taxonomy> --dbWGSDir <db_WGSdata_dir>

WGSValidate program is used to validate potential HGTs using WGS datasets.

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

";

#Read WGSValidation parameters
my $genome_file = "./kmerFilter/new.fasta";
my $HGT_id_file ="./screenHGT/HGT.id";
my $HGT_homologous_info_dir = "./screenHGT/tree/";
my $HGT_info_file = "./screenHGT/HGT.info";
my $type_file = "./splitDB/type.txt";
my $db_id_dir = "./splitDB/";
my $out_dir = "./WGSValidation/";
my $WGSflanking_length = 150;
my $depth = 10;
my $screen_distant = "kingdom";
my $screen_self = "phylum,class,order,species";
my $cdhit_threshold = 0.8;
my $mode = "Original";
my $db_dir;
my $db_info_file;
my $taxanomy_target;
my $db_WGSdata_dir;
my $help;
GetOptions(
    'genome=s'          => \$genome_file,
    'HGTId=s'		=> \$HGT_id_file,
    'HGTInfoDir=s'	=> \$HGT_homologous_info_dir,
    'HGTInfo=s'		=> \$HGT_info_file,
    'type=s'		=> \$type_file,
    'dbdir=s'           => \$db_dir,
    'dbInfo=s'		=> \$db_info_file,
    'dbIdDir=s'		=> \$db_id_dir,
    'taxo=s'            => \$taxanomy_target,
    'dbWGSDir=s'	=> \$db_WGSdata_dir,
    'outdir=s'		=> \$out_dir,
    'mode=s'		=> \$mode,
    'length=i'		=> \$WGSflanking_length,
    'depth=i'		=> \$depth,
    'help'      	=> \$help,
    'distant=s'		=> \$screen_distant,
    'self=s'            => \$screen_self,
    'idenCdHitEst=f'	=> \$cdhit_threshold
) or die $usage;

die $usage if defined($help);
die $usage if(!$db_dir);
die $usage if(!$db_info_file);
die $usage if(!$taxanomy_target);
die $usage if(!$db_WGSdata_dir);
die $usage if($mode ne "Original" && $mode ne "Strict");

my $flanking_length = $WGSflanking_length+350;
my $hit_file = "cov-hit.species.txt";
my $hit_file_WGS = "cov-hit.species.WGS.txt";
my ($genome_kingdom,$genome_phylum,$genome_class,$genome_order) = split(/,/,$taxanomy_target);
if($mode eq "Strict"){
    $hit_file = "cov-hit.species.strict.txt";
    $hit_file_WGS = "cov-hit.species.strict.WGS.txt";
    $HGT_id_file ="./screenHGT/modeStrict/HGT.id";
    $HGT_info_file ="./screenHGT/modeStrict/HGT.info";
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
unless($db_dir=~/\/$/){
    $db_dir.="/";
}
mkdir($out_dir);
my $out_dir_table = $out_dir."table/";
my $out_dir_faori = $out_dir."fa_ori/";
my $out_dir_fa = $out_dir."fa/";
my $out_dir_sam = $out_dir."sam/";
my $out_dir_cov = $out_dir."cov/";
my $out_dir_matrix = $out_dir."matrix/";
my $out_dir_afterWGS = $out_dir."afterWGS/";
mkdir($out_dir_table);
mkdir($out_dir_faori);
mkdir($out_dir_fa);
mkdir($out_dir_sam);
mkdir($out_dir_cov);
mkdir($out_dir_matrix);
mkdir($out_dir_afterWGS);

#Get homology sequences of candidate HGT sequences in other species
open(HGTID,$HGT_id_file)||die("screenHGT: error with opening $HGT_id_file\n");
while(<HGTID>){
    chomp();
    my $hgtid = $_;
    getHomoFA($hgtid,$HGT_homologous_info_dir,$db_dir,$flanking_length,$hit_file);
}
close HGTID;

#Extract HGT related species
extractSpecies($HGT_id_file,$HGT_homologous_info_dir,$db_info_file,$type_file,$out_dir_table,$hit_file);

#Filter HGT related species that require WGS validation
filterSpecies($HGT_info_file,$out_dir_table,$genome_kingdom);

#Acquire sequences that require WGS validation
acquireSeq($HGT_id_file,$HGT_homologous_info_dir,$out_dir_table,$out_dir_faori);
my $out_file_WGSspecies = $out_dir."table/WGS.species";
open(SPECIES,$out_file_WGSspecies) or die "WGSValidation.acquireSeq: error with opening $out_file_WGSspecies\n";
while(<SPECIES>){
    chomp();
    my $spe = $_;
    deleteDuplication($out_dir_faori,$out_dir_fa,$spe);
}
close SPECIES;

#Sequence alignemnt of homologous sequences and their WGSdata
open(SPECIES,$out_file_WGSspecies) or die "WGSValidation.acquireSeq: error with opening $out_file_WGSspecies\n";
while(<SPECIES>){
    chomp();
    my $spe = $_;
    seqAlignWGS($spe,$out_dir_fa,$out_dir_faori,$db_WGSdata_dir,$out_dir_sam,$out_dir_cov,$depth,$WGSflanking_length);
}
close SPECIES;
my $matrix = $out_dir_matrix."matrix.txt";
summaryCov($matrix,$out_dir_cov);

#Filter HGT homologous sequences that passed WGS validation
my $allid = $out_dir_matrix."all.id";
my $HGT_validated = $out_dir_matrix."validated.id";
my $HGT_notvalidated_ori = $out_dir_matrix."notvalidated.id.ori";
my $HGT_notvalidated_add = $out_dir_matrix."notvalidated.id.add";
my $HGT_notvalidated_mid = $out_dir_matrix."notvalidated.id.mid";
my $HGT_notvalidated = $out_dir_matrix."notvalidated.id";
my $comcat = "cat $out_dir_faori*id > $allid";
system($comcat);
my $WGSnoneedid = $out_dir_table."WGS_noneed.id";
validate($allid,$matrix,$HGT_validated,$HGT_notvalidated_ori);
del($allid,$HGT_notvalidated_add,$WGSflanking_length);
$comcat = "cat $HGT_notvalidated_ori $HGT_notvalidated_add > $HGT_notvalidated_mid";
system($comcat);
base::notin($WGSnoneedid,$HGT_notvalidated_mid,$HGT_notvalidated);
filterWGS($HGT_id_file,$HGT_homologous_info_dir,$db_info_file,$HGT_notvalidated,$hit_file,$hit_file_WGS);

#Filter potential HGTs using WGS datasets
#Calculate species distribution of candidate HGT sequences
screenHGT::calDistribution($out_dir_afterWGS,$screen_distant,$genome_order,$genome_class,$genome_phylum,$genome_kingdom,$HGT_id_file,$type_file,$HGT_homologous_info_dir,$hit_file_WGS);

#The minimum similarity of HGT in DRG is higher than the maximum similarity in CRG
my $distantid = $db_id_dir."e_".$screen_distant.".id";
screenHGT::distantIdenHclose($screen_distant,$screen_self,$hit_file_WGS,$db_id_dir,$type_file,$HGT_homologous_info_dir,$genome_file,$out_dir_afterWGS);

#Remove redundancy using cd-hit-est
my $dHcfa = $out_dir_afterWGS."distantHclose_".$screen_distant.".fa";
my $streeinfo = $out_dir_afterWGS."tree_".$screen_distant.".select.info";
screenHGT::rmRedundancy($dHcfa,$cdhit_threshold,$streeinfo,$out_dir_afterWGS);
}

sub getHomoFA{
my ($hgtid,$HGT_homologous_info_dir,$db_dir,$flanking_length,$hit_file) = @_;
my $file = $HGT_homologous_info_dir.$hgtid."/".$hit_file;
open(FILE,$file)||die("WGSValidation.getHomoFA: error with opening $file\n");
while(<FILE>){
    chomp();
    if(!($_ =~ /hit_species/)){
        my @arr = split(/\t/,$_);
        my $all = $_;
	my $flanking_st = 1;
        my $species = $HGT_homologous_info_dir.$hgtid."/hit/$arr[0].bed";
        my $species1 = $HGT_homologous_info_dir.$hgtid."/hit/$arr[0].bed.flanking";
        open(B, ">$species")||die("WGSValidation.getHomoFA: error with writing to $species\n");
        open(B1, ">$species1")||die("WGSValidation.getHomoFA: error with writing to $species1\n");
	if($arr[2]-$flanking_length > 1){
            $flanking_st = $arr[2]-$flanking_length;
        }
        my $flanking_en = $arr[3]+$flanking_length;
        print B "$arr[1]\t$arr[2]\t$arr[3]\t$arr[4]\n";
        my $out = $HGT_homologous_info_dir.$hgtid."/hit/$arr[0].fa";
        my $out1 = $HGT_homologous_info_dir.$hgtid."/hit/$arr[0].fa.flanking";
        open(OUT,">$out")||die("WGSValidation.getHomoFA: error with writing to $out\n");
        open(OUT1,">$out1")||die("WGSValidation.getHomoFA: error with writing to $out1\n");
        my $fa = $db_dir."$arr[0].fna";
        open(FA,$fa)||die("WGSValidation.getHomoFA: error with opening $fa\n");       
        my %hash = ();
        my %num = ();
        my $right = 0;
	my $faid = "";
        while(<FA>){
            chomp();
            if($_ =~ />([^\s]+)/){
                $faid=$1;
                if($faid eq $arr[1]){
                    $right=1;
                }elsif($right==0){
                    next;
                }else{
                    last;
                }
            }
            else{
                if($right){
                    $hash{$faid} = $hash{$faid}.$_;
                    $num{$faid} += length($_);
                }
            }  
        }
        my ($chr,$start,$end,$strand) = ($arr[1],$arr[2],$arr[3],$arr[4]);
        my ($re_start,$re_end) = (-$end,-$start);
        if($flanking_en>$num{$chr}){
            $flanking_en=$num{$chr};
        }
        print B1 "$arr[1]\t$flanking_st\t$flanking_en\t$arr[4]\n";
        my ($flanking_re_start,$flanking_re_end) = (-$flanking_en,-$flanking_st);
        my $length = $end-$start+1;
        my $flankinglength = $flanking_en-$flanking_st+1;
        my $HGT = "";
        my $flanking_HGT = "";
        $HGT = uc(substr($hash{$chr},$start-1,$end-$start+1));
        $flanking_HGT = uc(substr($hash{$chr},$flanking_st-1,$flankinglength));
        if($strand eq "-"){
            $HGT = uc(substr($hash{$chr},$re_start,$length));
            $HGT = reverse($HGT);
            $HGT =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;

            $flanking_HGT = uc(substr($hash{$chr},$flanking_re_start,$flankinglength));
            $flanking_HGT = reverse($flanking_HGT);
            $flanking_HGT =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
        }
        print OUT ">$chr|$start-$end|$strand\n$HGT\n";
        print OUT1 ">$chr|$flanking_st-$flanking_en|$strand\n$flanking_HGT\n";

        sub reverse_complement_IUPAC{         
        my ($dna) = @_;          
        my $revcomp = reverse($dna);
        $revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
        return $revcomp;
        }

    }
    close FA; close B; close B1; close OUT1; close OUT;
}
close FILE;

}

sub extractSpecies{
my ($HGT_id_file,$HGT_homologous_info_dir,$db_info_file,$db_type_file,$out_dir_table,$hit_file) = @_;
my $out_file_txt = $out_dir_table."speciesInfo.txt";
open(TXT,"> $out_file_txt") or die "WGSValidation.extractSpecies: error with writing to $out_file_txt\n";
open(NAME,$db_info_file) or die "WGSValidation.extractSpecies: error with opening $db_info_file\n";
my %name = ();
my @name = <NAME>;
foreach my $n (@name){
    chomp $n;
    my @n = split(/;|\t/,$n);
    $name{$n[0]} = $n[1];
}
open(TYPE,$db_type_file) or die "WGSValidation.extractSpecies: error with opening $db_type_file\n";
my %type = ();
my @type = <TYPE>;
foreach my $t (@type){
    chomp $t;
    my @t = split(/\t/,$t);
    my @gcfother = split(/_/,$t[0]);
    my $gcf = "$gcfother[0]_$gcfother[1]";
    $type{$gcf} = $t[1];
}
my %species = ();
foreach my $g (keys(%type)){
    $species{$g} = 0;
}
open(ID,$HGT_id_file) or die "WGSValidation.extractSpecies: error with opening $HGT_id_file\n";
my @id = <ID>;
foreach my $i (@id){
    chomp($i);
    my $single = $HGT_homologous_info_dir.${i}."/".$hit_file;
    open(SINGLE,$single) or die "WGSValidation.extractSpecies: error with opening $single\n";
    my @single = <SINGLE>;
    foreach my $s (@single){
        my @hit = split(/\t/,$s);
        my @gcfother = split(/_/,$hit[0]);
        my $gcf = "$gcfother[0]_$gcfother[1]";
        $species{$gcf} += 1;
        if($gcf ne "hit_species"){
            print TXT "$i\t$gcf\t$name{$gcf}\t$type{$gcf}\n";
        }
    }
}
close NAME; close TYPE; close ID; close TXT;
}

sub filterSpecies{
my ($HGT_info_file,$out_dir_table,$genome_kingdom) = @_;
my $out_file_txt = $out_dir_table."speciesInfo.txt";
my $out_file_WGSinfo = $out_dir_table."WGS.info";
my $out_file_WGSid = $out_dir_table."WGS.id";
my $out_file_WGSspecies = $out_dir_table."WGS.species";
my $out_file_WGSnoneedid = $out_dir_table."WGS_noneed.id";
open(TXT,$out_file_txt) or die "WGSValidation.extractSpecies: error with opening $out_file_txt\n";
open(INFO,$HGT_info_file) or die "WGSValidation.extractSpecies: error with opening $HGT_info_file\n";
open(WGSINFO,">$out_file_WGSinfo") or die "WGSValidation.extractSpecies: error with writing to $out_file_WGSinfo\n";
open(WGSID,">$out_file_WGSid") or die "WGSValidation.extractSpecies: error with writing to $out_file_WGSid\n";
open(WGSSPE,">$out_file_WGSspecies") or die "WGSValidation.extractSpecies: error with writing to $out_file_WGSspecies\n";
open(OUT,">$out_file_WGSnoneedid") or die "WGSValidation.extractSpecies: error with writing to $out_file_WGSnoneedid\n";

my @allkingdom = ("Bacteria","Viruses","Metazoa","Protozoa","Viridiplantae","Fungi");
foreach my $k (@allkingdom){
    if($k ne $genome_kingdom){
        push(@outkingdom,$k);
	%{$k} = ();
    }
}
my %id = ();
my %species = ();

while(<INFO>){
    chomp();
    @n = split(/\t/,$_);
    my $location = 5;
    foreach my $ok (@outkingdom){
        ${$ok}{$n[0]} = $n[$location];
	$location++;
    }
}
while(<TXT>){
    chomp();
    @a = split(/\t/,$_);
    if(exists($a[3]{$a[0]})){
        if($a[3]{$a[0]} < 10){
            $id{$a[1]}+=1;
            $species{$a[1]}=$a[2];
        }else{
            print OUT "$a[0]\t$a[2]\n";
        }
    }
}
foreach my $i (keys(%id)){
    print WGSINFO "$i\t$species{$i}\t$id{$i}\n";
    print WGSID "$i\n";
    print WGSSPE "$species{$i}\n";
}
close INFO; close TXT; close OUT; close WGSINFO; close WGSID; close WGSSPE;
}

sub acquireSeq{
my ($HGT_id_file,$HGT_homologous_info_dir,$out_dir_table,$out_dir_faori) = @_;
my $out_file_WGSid = $out_dir_table."WGS.id";
my $out_file_WGSspecies = $out_dir_table."WGS.species";
open(SPECIES,$out_file_WGSspecies) or die "WGSValidation.acquireSeq: error with opening $out_file_WGSspecies\n";
my @species = ();
while(<SPECIES>){
    chomp();
    push(@species,$_);
}
open(SID,$out_file_WGSid) or die "WGSValidation.acquireSeq: error with opening $out_file_WGSid\n";
my @sid = ();
while(<SID>){
    chomp();
    push(@sid,$_);
}
my $i = 0;
open(ORI,$HGT_id_file) or die "WGSValidation.acquireSeq: error with opening $HGT_id_file\n";
my %ori=();
while(<ORI>){
    chomp();
    $ori{$_}=0;
}
foreach $s (@species){
open(OUT, ">".$out_dir_faori.$s.".fa") or die "WGSValidation.acquireSeq: error with writing $out_dir_faori.$s\n";
open(OUT2, ">".$out_dir_faori.$s.".id") or die "WGSValidation.acquireSeq: error with writing $out_dir_faori.$s\n";
open(OUT3, ">".$out_dir_faori.$s.".fa.flanking") or die "WGSValidation.acquireSeq: error with writing $out_dir_faori.$s\n";
open(OUT4, ">".$out_dir_faori.$s.".id.flanking") or die "WGSValidation.acquireSeq: error with writing $out_dir_faori.$s\n";
my %hash = ();
opendir (DIR, $HGT_homologous_info_dir) or die "WGSValidation.acquireSeq: error with opening the directory $HGT_homologous_info_dir";
my @dir = readdir DIR;
foreach my $file (@dir){ 
    if(!($file =~ /hit/)){
    if(exists($ori{$file})){
    if(!exists($hash{$file})){
    $hash{$file} = 0;
    my $ndir = $HGT_homologous_info_dir.$file."/hit";
    opendir (NDIR, $ndir) or die "cannot open the directory";
    my @ndir = readdir NDIR;
    my %flank = ();
    foreach my $nfile (@ndir){
        if($nfile =~ /$sid[$i].*\.fa\.flanking/){
               open(REPEAT2, "< $ndir/$nfile") or die "WGSValidation.acquireSeq: error with opening $nfile";
               while(<REPEAT2>){
               chomp($_);
               if($_ =~ />(.*)/){
                  $flank{$sid[$i]}=$1;
               }
               }
        }
    }
    foreach my $nfile (@ndir){
        if($nfile =~ /$sid[$i].*\.fa/){
           if($nfile =~ /$sid[$i].*\.fa\.flanking/){
               open(REPEAT, "< $ndir/$nfile") or die "WGSValidation.acquireSeq: error with opening $nfile";
               while(<REPEAT>){
               chomp($_);
               print OUT3 "$_\n";
               if($_ =~ />(.*)/){
	           print OUT4 "$file\t$s\t$1\n";
               }
               }
           }elsif(!($nfile =~ /$sid[$i].*\.fa\.flanking/)){
               open(REPEAT1, "< $ndir/$nfile") or die "cannot open the $nfile";
               while(<REPEAT1>){
               chomp($_);
               print OUT "$_\n";
               if($_ =~ />(.*)/){
                  print OUT2 "$file\t$s\t$flank{$sid[$i]}\t$1\n";
               }
               }
           }
        }
    }
    }
    }
    }
}
$i = $i+1;
}
close SID; close SPECIES; close ORI; close DIR; close NDIR;
close REPEAT; close REPEAT2; close OUT; close OUT2; close OUT3; close OUT4;
}

sub deleteDuplication{
my($out_dir_faori,$out_dir_fa,$spe) = @_;
my $out_fa_ori = $out_dir_faori.$spe.".fa.flanking";
my $out_fa = $out_dir_fa.$spe.".fa";
open(FA,$out_fa_ori)||die("WGSValidation.deleteDuplication: error with opening $out_fa_ori\n");
open(OUT,">$out_fa")||die("WGSValidation.deleteDuplication: error with writing to $out_fa\n");
my %hash = ();
my $id = "";
while(<FA>){
    chomp();
    if($_ =~ />([^\s]+)/){
        $id = $1;
    }
    else{
        $hash{$id} = $_;
    }
}
my $h = "";
foreach $h (keys(%hash)){
    print OUT ">$h\n$hash{$h}\n";
}
close FA;close OUT;
}

sub seqAlignWGS{
use File::Basename;
use Cwd 'abs_path';
my ($spe,$out_dir_fa,$out_dir_faori,$db_WGS_dir,$out_dir_sam,$out_dir_cov,$depth,$WGSflanking_length) = @_;
my $out_fa = $out_dir_fa.$spe.".fa";
my $out_fa_index = $out_dir_fa.$spe;
my $combowtie2_build = "bowtie2-build $out_fa $out_fa_index";
system($combowtie2_build);

my @paired = ();
my @single = ();
my $db_WGS_dir_spe = $db_WGS_dir.$spe;
opendir(DIR, $db_WGS_dir) or die "WGSValidation.seqAlignWGS: error with opening the directory $db_WGS_dir";
my @dir = readdir DIR;
foreach my $ndir (@dir){
    if($ndir eq $spe){
	my %sra = ();
	opendir(NDIR, $db_WGS_dir_spe) or die "WGSValidation.seqAlignWGS: error with opening the directory $db_WGS_dir_spe";
	my @ndir = readdir NDIR;
	foreach my $file (@ndir){
	    if($file =~ /(.*)_1.fastq.gz$/ || $file =~ /(.*)_2.fastq.gz$/){
		if(not exists $sra{$1}){
		    @paired = (@paired,$1);
		    $sra{$1} = 1;
	        }
	    }elsif($file =~ /(.*).fastq.gz$/){
		if(not exists $sra{$1}){
		    @single = (@single,$1);
	       	    $sra{$1} = 1;
	        }
	    }
	}
    }
}

my $sam_dir = $out_dir_sam.$spe."/";
mkdir($sam_dir);
foreach my $sra (@paired){
    my $sam_dir_sra = $out_dir_sam.$spe."/".$sra.".sam";
    my $combowtie2 = "bowtie2 -x $out_fa_index -1 $db_WGS_dir_spe/${sra}_1.fastq.gz -2 $db_WGS_dir_spe/${sra}_2.fastq.gz -S $sam_dir_sra --no-unal --mm";
    system($combowtie2);
}
foreach my $sra (@single){
    my $sam_dir_sra = $out_dir_sam.$spe."/".$sra.".sam";
    my $combowtie2 = "bowtie2 -x $out_fa_index -U $db_WGS_dir_spe/${sra}.fastq.gz -S $sam_dir_sra --no-unal --mm";
    system($combowtie2);
}
my @allsra = (@paired,@single);
foreach my $sra (@allsra){
my $sam_dir_sra = $out_dir_sam.$spe."/".$sra;
my $comsamview = "samtools view -\@8 -b $sam_dir_sra.sam > $sam_dir_sra.bam";
my $comsamsort = "samtools sort $sam_dir_sra.bam -o $sam_dir_sra.sort";
my $comsamdepth = "samtools depth -a $sam_dir_sra.sort > $sam_dir_sra.txt";
my $comrm = "rm -rf $sam_dir_sra.sam $sam_dir_sra.bam $sam_dir_sra.sort";
system($combowtie2);
system($comsamview);
system($comsamsort);
system($comsamdepth);
system($comrm);
}

my $cov_dir = $out_dir_cov.$spe."/";
mkdir($cov_dir);
my $faid = $out_dir_faori.$spe.".id";
foreach my $sra (@allsra){
my $samdepth = $sam_dir.$sra.".txt";
my $cov = $cov_dir.$sra.".cov";
my $covnormal = $cov_dir.$sra.".cov.normal";
my @args = ($samdepth,$faid,$cov,$covnormal,$depth,$WGSflanking_length);
my $module_dir = dirname(abs_path(__FILE__));
system("python","$module_dir/ttest.py",@args) == 0 or die "WGSValidation.seqAlignWGS: system call python failed: $?";
}
}

sub summaryCov{
my ($matrix,$out_dir_cov) = @_;
open(MATRIX,">$matrix")||die("WGSValidation.summaryCov: error with writing to $matrix\n");
opendir(DIR,$out_dir_cov) or die "WGSValidation.summaryCov: error with opening the directory $out_dir_cov";
my @dir = readdir DIR;
foreach my $species (@dir){
    opendir(NDIR,$out_dir_cov.$species) or die "WGSValidation.summaryCov: error with opening the directory $species";
    my @files = readdir NDIR;
    foreach my $samples (@files){
	if($samples =~ /([^\s]+)\.cov$/){
	    my $sample = $1;
	    my $file = $out_dir_cov.$species."/".$samples;
	    open(FILE,$file) or die "WGSValidation.summaryCov: error with opening $file";
	    while(<FILE>){
                chomp($_);
        	my @data = split(/\s+/,$_);
        	print MATRIX "$species\t$data[0]\t$sample\t$data[1]\n";
    	    }
        }
    }
}
}

sub validate{
my ($allid,$matrix,$HGT_validated,$HGT_notvalidated_ori) = @_;
open(ALL,$allid) or die "WGSValidation.validate: error with opening $allid\n";
open(VAL,$matrix) or die "WGSValidation.validate: error with opening $matrix\n";
open(OUT,">$HGT_validated") or die "WGSValidation.validate: error with writing to $HGT_validated\n";
open(OUT1,">$HGT_notvalidated_ori") or die "WGSValidation.validate: error with writing to $HGT_validated_ori\n";
my %aval=();
my $avaled=();
my %val=();
my $valed=();

while(<VAL>){
    chomp();
    my @vdata=split(/\t/,$_);
    my $vkey="$vdata[0]\t$vdata[1]";
    if(!exists($val{$vkey})){
        $val{$vkey}=0;
        $valed{$vkey}=0;
    }
    $val{$vkey}+=1;
    if($vdata[3]>=0.8){
        $valed{$vkey}+=1;
    }
}        
while(<ALL>){
    chomp();
    my @data=split(/\t/,$_);
    my $key="$data[1]\t$data[3]";
    my $akey=$data[0];
    $aval{$akey}+=$val{$key};
    $avaled{$akey}+=$valed{$key};
    if($valed{$key}<2){
        print OUT1 "$akey\t$data[1]\n";
    }
}
foreach my $k (keys(%aval)){
    print OUT "$k\t$avaled{$k}\t$aval{$k}\n";
}
close ALL; close VAL; close OUT1; close OUT;
}

sub del{
my ($allid,$HGT_notvalidated_add,$WGSflanking_length) = @_;
open(OUT,">$HGT_notvalidated_add") or die "WGSValidation.del: error with writing to $HGT_notvalidated_add\n";
open(ID,$allid) or die "WGSValidation.del: error with opening $allid\n";
while(<ID>){
    chomp();
    my @arr=split(/\t/,$_);
    my @flanking=split(/-|\|/,$arr[2]);
    my @ori=split(/-|\|/,$arr[3]);
    if($flanking[1]>$ori[1]-$WGSflanking_length || $flanking[2]<$ori[2]+$WGSflanking_length){
        print OUT "$arr[0]\t$arr[1]\n";
    }
}
close OUT; close ID;
}

sub filterWGS{
my ($HGT_id_file,$HGT_homologous_info_dir,$db_info_file,$HGT_notvalidated,$hit_file,$hit_file_WGS) = @_;
open(ID2SPECIES,$db_info_file)||die("WGSValidation.filterWGS: error with opening $db_info_file\n");
my %id = ();
while(<ID2SPECIES>){
    chomp();
    my @arr = split(/;/,$_);
    $id{$arr[1]} = $arr[0];
}
open(DELETE,$HGT_notvalidated)||die("WGSValidation.filterWGS: error with opening $HGT_notvalidated\n");
my %delete=();
while(<DELETE>){
    chomp();
    my @arr = split(/\t/,$_);
    $delete{$arr[0]}.="\t$id{$arr[1]}";
}
open(ORI,$HGT_id_file) or die "WGSValidation.filterWGS: error with opening $HGT_id_file\n";
my @id=<ORI>;
foreach $i (@id){
    chomp($i);
    my $file = $HGT_homologous_info_dir.$i."/".$hit_file;
    my $out = $HGT_homologous_info_dir.$i."/".$hit_file_WGS;
    open(OUT,">$out")||die("WGSValidation.filterWGS: error with writing to $out\n");
    open(F,$file)||die("WGSValidation.filterWGS: error with opening $file\n");
    while(<F>){
        chomp();
        my @arr = split(/_/,$_);
        if(!(exists($delete{$i}) && $delete{$i}=~ /GCF_$arr[1]/)){
            print OUT "$_\n";
        }
    }
    close F; close OUT;
}
close ORI; close DELETE; close ID2SPECIES;
}

1;
