#!/usr/bin/perl
package kmerFilter;

sub kmerFilter{
use strict;
use warnings;
use Getopt::Long;

my $usage="\nUsage: GIFEHGT kmerFilter [options] --fasta <genome_file>

kmerFilter program is used to select most different fragments with target genome by calculate the distance of kmer frequecies.

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
     
    --kmerPer	<int>		The percentage of fragments selected that are most different from the target genome. (100 in standard mode and 20 in accelerated mode)

    --trfPer	<int>           The percentage of sequences overlapping with simple repeats, and sequences above this percentage are judged as repeats and deleted as a whole. (50)
     
    --trfPara	<string>	The parameters of program trf (a comma-separated string), including matching weight,mismatching penaltyindel penalty,match probability,indel probability,minimum alignment score to report,maximum period size to report,maximum TR length expected. (2,5,7,80,10,50,1000,3)

";

#read kmerFilter parameters
my $out_dir = "./kmerFilter/";
my $mode = "standard";
my $split_len=1000;
my $split_step=800;
my $kmerFilter_kmer=4;
my $kmerFilter_per=100;
my $remove_per = 0.5;
my ($remove_match,$remove_Mismatch,$remove_Delta,$remove_PM,$remove_PI,$remove_Minscore,$remove_MaxPeriod,$remove_len) = (2,5,7,80,10,50,1000,3);
my $genome_file;
my $kmerFilter_per_input=20;
my $help;
my $remove_parameter;

GetOptions(
    'fasta=s'	=> \$genome_file,
    'mode=s'    => \$mode,
    'outdir=s'  => \$out_dir,
    'help'      => \$help,
    'length=i'  => \$split_len,
    'step=i'    => \$split_step,
    'kmer=i'    => \$kmerFilter_kmer,
    'kmerPer=i' => \$kmerFilter_per_input,
    'trfPer=i'  => \$remove_per,
    'trfPara=s' => \$remove_parameter
) or die $usage;
($remove_match,$remove_Mismatch,$remove_Delta,$remove_PM,$remove_PI,$remove_Minscore,$remove_MaxPeriod,$remove_len) = split(/,/, $remove_parameter) if $remove_parameter;

die $usage if defined($help);
die $usage if(!$genome_file);

#Check existence of output directory
if(-e $out_dir){
    print("Warning: output directory \"$out_dir\" already exists. Existing files will be overwritten.\n");
    `rm -rf $out_dir`;
}

#Adjust kmerFilter_per
if($mode eq "accelerated"){
    $kmerFilter_per=$kmerFilter_per_input;
}

#Adjust directory names and create output directory
unless($out_dir=~/\/$/){
    $out_dir.="/";
}
mkdir($out_dir);
my $new_genome_file=$out_dir."new.fasta";
my $split_fragments_file=$out_dir."split.fasta";
my $kmer_genome_file=$out_dir."kmer.txt";
my $kmer_fragments_file=$out_dir."fragments_kmer.txt";
my $kmer_distance_file=$out_dir."kmer_distance.txt";
my $kmer_filter_file=$out_dir."kmerFilter.fasta";
my $kmer_filter_name="kmerFilter.fasta";
my $trf_out_name = "kmerFilter";
my $trf_out_file = $out_dir.$trf_out_name.".trf.mask";
my $outbed = $out_dir.$trf_out_name."_removeSimRep.bed";
my $outfa = $out_dir.$trf_out_name."_removeSimRep.fa";

my @keys = ();
@keys = generateKmer($kmerFilter_kmer,@keys);

#Convert fasta files from a multi-line format to a single-line format
convert($genome_file,$new_genome_file);

#Split the target genome into fragments
splittofragments($new_genome_file,$split_len,$split_step,$split_fragments_file);

#Calculate kmer frequencies of target genome
kmerGenome($new_genome_file,$kmerFilter_kmer,$kmer_genome_file,\@keys);

#Calculate kmer frequencies of split fragments
kmerFragments($split_fragments_file,$kmerFilter_kmer,$kmer_fragments_file,\@keys);

#Calculate distance between kmer frequencies of split fragments and that of target genome
calculateDistance($kmer_genome_file,$kmer_fragments_file,$kmer_distance_file);

#Select most different fragments with target genome
selectFragments($new_genome_file,$kmer_distance_file,$kmerFilter_per,$kmer_filter_file);

#Remove simple repeat sequences of input file using trf
#Run trf to mask the simple repeat sequences
runtrf($kmer_filter_file,$kmer_filter_name,$trf_out_file,$remove_match,$remove_Mismatch,$remove_Delta,$remove_PM,$remove_PI,$remove_Minscore,$remove_MaxPeriod,$remove_len);
#Remove fragments with simple repeat sequences
removeSim($trf_out_file,$remove_per,$outbed,$new_genome_file,$outfa);

}

#Convert fasta files from a multi-line format to a single-line format
sub convert{
my ($genome_file,$mid_file)=@_;

open(FA,$genome_file)||die("splitGenome.convert: error with opening the fasta file of target genome\n");
open(OUT,">$mid_file")||die("splitGenome.convert: error with writing to the output directory\n");

my $current_id = '';
my $current_sequence = '';

while (my $line = <FA>) {
    chomp $line;
    if ($line =~ /^>(.+)/) {
        if ($current_id) {
            print OUT ">$current_id\n$current_sequence\n";
        }
        $current_id = $1;
        $current_sequence = '';
    } else {
        $current_sequence .= $line;
    }
}

#for last sequence
if ($current_id) {
    print OUT ">$current_id\n$current_sequence\n";
}

close FA;close OUT;
}

#Split the target genome into fragments
sub splittofragments{
my ($mid_file,$split_len,$split_step,$out_file)=@_;

# open file
open(FA,$mid_file)||die("splitGenome.splittofragments: error with opening the fasta file of target genome\n");
open(OUT,">$out_file")||die("splitGenome.splittofragments: error with writing to the output directory\n");

# split the target genome to fragments
my $id = "";
while(<FA>){
    chomp();
    if($_ =~ />([^\s]+)/){
        $id = $1;
    }
    else{
        my $length = length($_);
        for(my $i=0;$i<$length-$split_len;){
            my $seq = substr($_,$i,$split_len);
            my ($start,$end) = ($i+1,$i+$split_len);
            print OUT ">$id-$start-$end\n$seq\n";
            $i += $split_step
        }
    }
}

close FA;close OUT;
}

#Calculate kmer frequencies of target genome
sub kmerGenome{
my ($genome_file,$kmerFilter_kmer,$kmer_genome_file,$keys)=@_;

if($genome_file =~ /\.gz/){
    open(GENOME,"gzip -dc $genome_file|")||die("kmerFilter.kmerGenome: error with opening the fasta file of target genome\n");
}
else{
    open(GENOME,$genome_file)||die("kmerFilter.kmerGenome: error with opening the fasta file of target genome\n");
}
open(OUT,">$kmer_genome_file")||die("kmerFilter.kmerGenome: error with writing to the output directory\n");

my %hash = ();
my @keys = @{$keys};

foreach my $key(@keys){
    $hash{$key} = 0;
}

my $sum = 0;
while(<GENOME>){
    chomp();
    if($_ !~ />/){
        my $seq = $_;
        my $len = length($seq);
        for(my $i=0;$i<$len-$kmerFilter_kmer+1;$i++){
            my $base = uc(substr($seq,$i,$kmerFilter_kmer));
            if(exists($hash{$base})){
                $hash{$base} += 1;
                $sum += 1;
            }
        }
    }
}

print OUT "wholeGenome\t";
foreach my $key(@keys){
    $hash{$key} = sprintf("%0.4f",$hash{$key}/$sum);
    print OUT "$hash{$key}\t";
}
print OUT "\n";

close GENOME; close OUT;
}

#Calculate kmer frequencies of split fragments
sub kmerFragments{
my ($fragments_file,$kmerFilter_kmer,$kmer_fragments_file,$keys)=@_;

open(Fragment,$fragments_file)||die("kmerFilter.kmerFragments: error with opening the file of split fragments\n");
open(OUT,">$kmer_fragments_file")||die("kmerFilter.kmerFragments: error with writing to the output directory\n");

my %hash = ();
my @keys = @{$keys};
foreach my $key(@keys){
    $hash{$key} = 0;
}

my $sum = 0;
my $id = "";
while(<Fragment>){
    chomp();
    if($_ =~ />([^\s]+)/){
        $id = $1;
    }
    if($_ !~ />/){
        my $seq = $_;
        my $len = length($seq);
        for(my $i=0;$i<$len-$kmerFilter_kmer+1;$i++){
            my $base = uc(substr($seq,$i,$kmerFilter_kmer));
            if(exists($hash{$base})){
                $hash{$base} += 1;
                $sum += 1;
            }
        }
        if($sum > 0){
            print OUT "$id\t";
            foreach my $key(@keys){
                $hash{$key} = sprintf("%0.4f",$hash{$key}/$sum);
                print OUT "$hash{$key}\t";
            }
            print OUT "\n";
        }
        foreach my $key(@keys){
            $hash{$key} = 0;
        }
        $sum = 0;
    }
}
close Fragment; close OUT;
}

#Calculate distance between kmer frequencies of split fragments and that of target genome
sub calculateDistance{
my ($kmer_genome_file,$kmer_fragments_file,$kmer_distance_file)=@_;

open(BACKGROUND,$kmer_genome_file)||die("kmerFilter.calculateDistance: error with opening the file of target genome\n");
open(TARGET,$kmer_fragments_file)||die("kmerFilter.calculateDistance: error with opening the file of split fragments\n");
open(OUT,">$kmer_distance_file")||die("kmerFilter.calculateDistance: error with writing to the output directory\n");

my @genome = ();
while(<BACKGROUND>){
    chomp();
    @genome = split(/\s+/,$_);
    last;
}

print OUT "region\tdistance\n";
while(<TARGET>){
    chomp();
    my @data = split(/\s+/,$_);
    my $distance = 0;
    my $size = @genome;
    for(my $i=1;$i<$size;$i++){
        $distance += ($data[$i]-$genome[$i])**2;
    }
    $distance = sqrt($distance);
    print OUT "$data[0]\t$distance\n";
}

close BACKGROUND;close TARGET;close OUT;
}

#Select most different fragments with target genome
sub selectFragments{
use lib '.';
use base;
my ($genome_file,$kmer_distance_file,$kmerFilter_per,$out_file)=@_;

open(IN,$kmer_distance_file)||die("kmerFilter.selectFragments: error with opening the distance file\n");

my $totalnum = `grep -v region $kmer_distance_file | wc -l`;
chomp($totalnum);
my $midnum = `expr 100 / $kmerFilter_per `;
chomp($midnum);
my $selectnum = `expr $totalnum / $midnum`;
chomp($selectnum);
my $outputinfo = `grep -v region $kmer_distance_file |sort -rnk 2 | head -$selectnum |awk '{print $1}'`;
$outputinfo =~ s/-/\t/g;
getSeq($genome_file,$outputinfo,$out_file);

close IN; 
}

sub generateKmer{
    my ($order,@previous) = @_;
    my @append = ();
    
    if(@previous == 0){
        @append = ('A','T','C','G');
    }
    else{
        my @base = ('A','T','C','G');
        foreach my $item(@previous){
            foreach my $char(@base){
                my $item2 = $item.$char;
                push(@append,$item2);
            }
        }
    }
    
    $order -= 1;
    if($order == 0){
        return @append;
    }
    else{
        generateKmer($order,@append);
    }
}

sub getSeq{
my ($fa,$bed,$out)= @_;
open(FA,$fa)||die("error with fasta\n");
open(OUT,">$out")||die("error with out\n");

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

my @bed = split(/\n/, $bed);
foreach $b (@bed){
    chomp($b);
    my @arr = split(/\s+/,$b);
    my ($chr,$start,$end) = ($arr[0],$arr[1],$arr[2]);
    my $HGT = uc(substr($hash{$chr},$start-1,$end-$start+1));
    print OUT ">$chr-$start-$end\n$HGT\n";
}

close FA;close OUT;
}

#Run trf to mask the simple repeat sequences
sub runtrf{
my ($fasta_file,$fasta_name,$trf_out_file,$remove_match,$remove_Mismatch,$remove_Delta,$remove_PM,$remove_PI,$remove_Minscore,$remove_MaxPeriod,$remove_len)=@_;

open(FA,$fasta_file)||die("removeSimRep.runtrf: error with opening the fasta file\n");

my $comtrf = "trf $fasta_file $remove_match $remove_Mismatch $remove_Delta $remove_PM $remove_PI $remove_Minscore $remove_MaxPeriod -l $remove_len -m -h";
system($comtrf);
my $commv = "mv $fasta_name.$remove_match.$remove_Mismatch.$remove_Delta.$remove_PM.$remove_PI.$remove_Minscore.$remove_MaxPeriod.mask $trf_out_file";
my $comrm = "rm -rf $fasta_name.$remove_match.$remove_Mismatch.$remove_Delta.$remove_PM.$remove_PI.$remove_Minscore.$remove_MaxPeriod.dat";
system($commv);
system($comrm);

close FA;
}

sub removeSim{
my ($trf_out_file,$remove_per,$outbed,$genome_file,$outfa)=@_;

trfFilter($trf_out_file,$remove_per,$outbed);
getSeqfile($genome_file,$outbed,$outfa);
}

#Remove fragments with simple repeat sequences
sub trfFilter{
my ($trf_out_file,$remove_per,$outbed)=@_;

open(MASK,$trf_out_file)||die("removeSimRep.removeSim.trfFilter: error with opening the trf outfile\n");
open(OUT,">$outbed")||die("removeSimRep.removeSim.trfFilter: error with opening the output file\n");

my %seq=();

my $s="";
while(<MASK>){
    chomp();
    if($_ =~ />(?<chr>.*)-(?<start>.*)-(?<end>\d+)/){
        $s="$+{chr}-$+{start}-$+{end}";
        $seq{$s}="";
    }else{
        $seq{$s}=$seq{$s}.$_;
    }
    #print "seq\t$s\t$seq{$s}\n";
}

foreach my $k (keys %seq){
    @split=split(/-/,$k);
    my $length=$split[2]-$split[1]+1;
    my $repeat=($seq{$k}=~s/N/N/g);
    my $cov=$repeat/$length;
    #print "$split[0]-$split[1]-$split[2]\t$seq{$k}\n$repeat\t$length\t$cov\t$remove_per\n";
    if($cov < $remove_per){
        print OUT "$split[0]\t$split[1]\t$split[2]\n";
    }
}

close MASK; close OUT;
}

sub getSeqfile{
my ($fa,$bed,$out)= @_;
open(FA,$fa)||die("getSeq: error with fasta $fa\n");
open(BED,$bed)||die("getSeq: error with bed $bed\n");
open(OUT,">$out")||die("getSeq: error with out $out\n");

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

while(<BED>){
    chomp();
    my @arr = split(/\s+/,$_);
    my ($chr,$start,$end) = ($arr[0],$arr[1],$arr[2]);
    my $HGT = uc(substr($hash{$chr},$start-1,$end-$start+1));
    print OUT ">$chr-$start-$end\n$HGT\n";
}

close FA;close BED;close OUT;
}

1;
