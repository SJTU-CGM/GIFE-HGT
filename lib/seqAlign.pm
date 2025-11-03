#!/usr/bin/perl
package seqAlign;

sub seqAlign{
use strict;
use warnings;
use Getopt::Long;

my $usage="\nUsage: GIFEHGT seqAlign [options] --db <db_genome_dir>

SeqAlignprogram is used to sequence alignment using LASTZ.

Necessary input description:

  db_genome_dir		<string>	The directory of all genomes saved.

Options (defaults in parentheses):

    --help				Print this usage page.
  
  Input:

    --fragment		<string>	The fasta file of filtered fragments of target genome. (./kmerFilter/kmerFilter_removeSimRep.fa)

    --genome		<string>	The fasta file of target genome. (./kmerFilter/new.fasta)

    --dbinfodir		<string>	The directory that all genome information saved. (./splitDB/)

  Output:

    --outdir		<string>	The result files will be output to this directory.

  Other parameters:

    --suffix		<string>	The suffix of genome file in database. (fna)

    --distant		<string>        The distantly related group. kingdom or phylum can be chosen. (kingdom)

    --identity		<float>		The threshold of identity. Two sequences above this threshold are considered to be similar. (0.5)

    ";

#read seqAlign parameters
my $fragments_file = "./kmerFilter/kmerFilter_removeSimRep.fa";
my $genome_file = "./kmerFilter/new.fasta";
my $db_info_dir = "./splitDB/";
my $out_dir = "./seqAlign/";
my $genome_suffix = "fna";
my $iden_threshold = 0.5;
my $screen_distant = "kingdom";
my $db_genome_dir;
my $help;
GetOptions(
    'fragment=s'	=> \$fragments_file,
    'genome=s'		=> \$genome_file,
    'db=s'		=> \$db_genome_dir,
    'dbinfo=s'		=> \$db_info_dir,
    'outdir=s' 		=> \$out_dir,
    'help'      	=> \$help,
    'suffix=s'  	=> \$genome_suffix,
    'distant=s'		=> \$screen_distant,
    'identity=f'	=> \$iden_threshold
) or die $usage;

die $usage if defined($help);
die $usage if (!$db_genome_dir);

#Check existence of output directory
if(-e $out_dir){
    print("Warning: output directory \"$out_dir\" already exists. Existing files will be overwritten.\n");
    `rm -rf $out_dir`;
}

#Adjust directory names and create output directory
unless($out_dir=~/\/$/){
    $out_dir.="/";
}
unless($db_genome_dir=~/\/$/){
    $db_genome_dir.="/";
}
unless($db_info_dir=~/\/$/){
    $db_info_dir.="/";
}
my $out_dir_distant = $out_dir."distant/";
my $out_dir_close = $out_dir."close/";
my $out_dir_close_fa = $out_dir_close."fa/";
mkdir($out_dir);
mkdir($out_dir_distant);
mkdir($out_dir_close);
mkdir($out_dir_close_fa);

my $closeid = $db_info_dir.$screen_distant.".id";
my $distantid = $db_info_dir."e_".$screen_distant.".id";

#Sequence alignment with species in DRG
sequenceAlignment($distantid,$out_dir_distant,$fragments_file,$db_genome_dir,$genome_suffix);

#Acquire target genome sequences which can align to species in DRG
acquireSeqDRG($out_dir_distant,$out_dir_close_fa,$genome_file);

#Sequence alignment with species in CRG
my $genome_aligndistant_file = $out_dir_close_fa."distant.fa";
sequenceAlignment($closeid,$out_dir_close,$genome_aligndistant_file,$db_genome_dir,$genome_suffix);

#Move file in DRG and CRG to the sam directory
moveFile($out_dir);

}

#Sequence alignment and convert format
sub sequenceAlignment{
my ($genomeid,$out_dir,$genome_file,$db_genome_dir,$genome_suffix) = @_;
my $out_dir_axt = $out_dir."axt/";
my $out_dir_cov = $out_dir."cov/";
my $out_dir_iden = $out_dir."iden/";
mkdir($out_dir_axt);
mkdir($out_dir_cov);
mkdir($out_dir_iden);
open(ID,$genomeid)||die("seqAlign.sequenceAlignment: error with opening $genomeid\n");
while(<ID>){
    chomp();
    my $id=$_;
    my $axtfile = $out_dir_axt.$id.".axt";
    my $covfile = $out_dir_cov.$id.".cov";
    my $bedfile = $out_dir_cov.$id.".bed";
    my $mergefile = $out_dir_cov.$id.".bed.merge";
    my $idenfile = $out_dir_iden.$id.".bed";
    #Sequence alignment using LASTZ
    lastz($genome_file,$db_genome_dir,$id,$genome_suffix,$axtfile);
    #convert axt format to cov format
    axt2cov($axtfile,$covfile);
    #convert cov format to bed format
    cov2bed($covfile,$bedfile);
    #merge bed
    mergeBed($bedfile,$mergefile);
    #select frangments with identity more than threshold
    selectIden($mergefile,$idenfile,$iden_threshold);
}
}

#Acquire target genome sequences which can align to species in DRG
sub acquireSeqDRG{
my ($out_dir_distant,$out_dir_close_fa,$genome_file) = @_;
my $aligndistantabed = $out_dir_distant."iden/*bed";
my $aligndistant_bed = $out_dir_close_fa."distant.bed";
my $aligndistant_bed_sort = $out_dir_close_fa."distant.bed.sort";
my $genome_aligndistant_file = $out_dir_close_fa."distant.fa";
my $comcat = "cat $aligndistantabed > $aligndistant_bed";
system($comcat);
base::sort($aligndistant_bed,$aligndistant_bed_sort);
base::merge($aligndistant_bed_sort,$aligndistant_bed);
my $comrm = "rm -rf $aligndistant_bed_sort";
system($comrm);
base::getSeq($genome_file,$aligndistant_bed,$genome_aligndistant_file);
}

#Move file in DRG and CRG to the sam directory
sub moveFile{
my ($out_dir) = @_;
my $out_dir_all = $out_dir."all/";
my $out_dir_all_iden = $out_dir."all/iden/";
my $out_dir_all_cov = $out_dir."all/cov/";
mkdir($out_dir_all);
mkdir($out_dir_all_iden);
mkdir($out_dir_all_cov);
my $out_dir_distant_iden = $out_dir."distant/iden/*";
my $out_dir_close_iden = $out_dir."close/iden/*";
my $commvidenc = "mv $out_dir_close_iden $out_dir_all_iden";
my $commvidend = "mv $out_dir_distant_iden $out_dir_all_iden";
system($commvidenc);
system($commvidend);
my $out_dir_distant_cov = $out_dir."distant/cov/*merge";
my $out_dir_close_cov = $out_dir."close/cov/*merge";
my $commvcovc = "mv $out_dir_close_cov $out_dir_all_cov";
my $commvcovd = "mv $out_dir_distant_cov $out_dir_all_cov";
system($commvcovc);
system($commvcovd);
}

sub lastz{
my ($genome_file,$db_dir,$id,$genome_suffix,$axtfile) = @_;
my $input1 = $genome_file."[multiple]";
my $input2 = $db_dir.${id}.".$genome_suffix";
my $comlastz = "lastz_32 $input1 $input2 --format=axt+ --output=$axtfile --ambiguous=iupac";
system ($comlastz);
}

sub axt2cov{
my ($axtfile,$covfile) = @_;
open(IN,$axtfile)||die("seqAlign.axt2cov: error with opening $axtfile\n");
open(OUT, ">$covfile")||die("seqAlign.axt2cov: error with writing to $covfile\n");
    
my ($chr1,$start1,$end1,$chr2,$start2,$end2,$strand2,$score,$query,$target,$identity) = ("","","","","","","","","","","");
my $index = 0; ##line index
while(<IN>){
    chomp();
    if($_ !~ /#/ && $_ =~ /[^\s]/){
        if($index%3 == 0){
            my @info = split(/\s+/);
	    ($chr1,$start1,$end1,$chr2,$start2,$end2,$strand2,$score) = ($info[1],$info[2],$info[3],$info[4],$info[5],$info[6],$info[7],$info[8]);
        }elsif($index%3 == 1){
            $query = uc($_);
        }else{
            $target = uc($_);
            my $len1 = length($query);
	    my $len2 = length($target);
            my $score = 0;
            for(my $i=0;$i<$len1 && $i<$len2;$i++){
                my $chr1 = substr($query,$i,1);
		my $chr2 = substr($target,$i,1);
                if($chr1 eq $chr2){
                    $score += 1;
                }
            }
            $identity = sprintf("%0.3f",$score/$len1);
            print OUT "$chr1\t$start1\t$end1\t$chr2\t$start2\t$end2\t$strand2\t$score\t$identity\n";
        }
        $index += 1;
    }
}

close IN; close OUT;

my $covfile_sort = $covfile.".sort";
`sort -k1,1 -k2,2n $covfile > $covfile_sort`;
`mv $covfile_sort $covfile`;
}

sub cov2bed{
my ($covfile,$bedfile) = @_;
open(IN,$covfile)||die("seqAlign.cov2bed: error with opening $covfile\n");
open(OUT, ">$bedfile")||die("seqAlign.cov2bed: error with writing to $bedfile\n");
my ($chr,$start,$end,$chr_0,$start_0,$end_0,$strand,$identity,$length) = ();
while(<IN>){
    chomp($_);
    my @arr = split(/\s+/,$_);
    my $chr1 = $arr[0];
    my $start1 = $arr[1];
    my $end1 = $arr[2];
    my $identity1 = $arr[8];
    my $length1 = $end1-$start1+1;
    my $chr1_0 = $arr[3];
    my $start1_0 = $arr[4];
    my $end1_0 = $arr[5];
    my $strand1 = $arr[6];
    if($chr eq ""){
        ($chr,$start,$end,$chr_0,$start_0,$end_0,$strand,$identity,$length) = ($chr1,$start1,$end1,$chr1_0,$start1_0,$end1_0,$strand1,$identity1,$length1);
    }
    elsif($chr1 ne $chr){
        @data = split(/-/,$chr);
        $nstart=$data[1]+$start-1;
        $nend=$data[1]+$end-1;
        print OUT "$data[0]\t$nstart\t$nend\t$chr_0\t$start_0\t$end_0\t$strand\t$identity\n";
        ($chr,$start,$end,$chr_0,$start_0,$end_0,$strand,$identity,$length) = ($chr1,$start1,$end1,$chr1_0,$start1_0,$end1_0,$strand1,$identity1,$length1);
    }
    else{
        if($start1 > $end + 1){
            @data = split(/-/,$chr);
            $nstart=$data[1]+$start-1;
            $nend=$data[1]+$end-1;
            print OUT "$data[0]\t$nstart\t$nend\t$chr_0\t$start_0\t$end_0\t$strand\t$identity\n";
            ($start,$end,$chr_0,$start_0,$end_0,$strand,$identity,$length) = ($start1,$end1,$chr1_0,$start1_0,$end1_0,$strand1,$identity1,$length1);
        }else{
            if(($length1*$identity1 > $length*$identity && $identity1 > 0.95*$identity) || $identity < 0.95*$identity1){
                ($start,$end,$chr_0,$start_0,$end_0,$strand,$identity,$length) = ($start1,$end1,$chr1_0,$start1_0,$end1_0,$strand1,$identity1,$length1);
            }
        }
    }
}
@data = split(/-/,$chr);
$nstart=$data[1]+$start-1;
$nend=$data[1]+$end-1;
print OUT "$data[0]\t$nstart\t$nend\t$chr_0\t$start_0\t$end_0\t$strand\t$identity\n";

close IN; close OUT;

my $bedfile_sort = $bedfile.".sort";
`sort -k1,1 -k2,2n $bedfile > $bedfile_sort`;
`mv $bedfile_sort $bedfile`;
}

sub mergeBed{
my ($bedfile,$mergefile) = @_;
open(IN,$bedfile)||die("seqAlign.mergeBed: error with opening $bedfile\n");
open(OUT, ">$mergefile")||die("seqAlign.mergeBed: error with writing to $mergefile\n");

my ($chr,$start,$end,$chr_0,$start_0,$end_0,$strand,$identity,$length) = ("",-1,-1,"",0,0,"",0,0);
while(<IN>){
    chomp($_);
    my @arr = split(/\t/,$_);
    my $chr1 = $arr[0];
    my $start1 = $arr[1];
    my $end1 = $arr[2];
    my $identity1 = $arr[7];
    my $length1 = $end1-$start1+1;
    my $chr1_0 = $arr[3];
    my $start1_0 = $arr[4];
    my $end1_0 = $arr[5];
    my $strand1 = $arr[6];
    if($chr eq ""){
        ($chr,$start,$end,$chr_0,$start_0,$end_0,$strand,$identity,$length) = ($chr1,$start1,$end1,$chr1_0,$start1_0,$end1_0,$strand1,$identity1,$length1);
    }elsif($chr1 ne $chr){
        print OUT "$chr\t$start\t$end\t$chr_0\t$start_0\t$end_0\t$strand\t$identity\n";
        ($chr,$start,$end,$chr_0,$start_0,$end_0,$strand,$identity,$length) = ($chr1,$start1,$end1,$chr1_0,$start1_0,$end1_0,$strand1,$identity1,$length1);
    }
    else{
        if($start1 > $end + 1){
            print OUT "$chr\t$start\t$end\t$chr_0\t$start_0\t$end_0\t$strand\t$identity\n";
            ($start,$end,$chr_0,$start_0,$end_0,$strand,$identity,$length) = ($start1,$end1,$chr1_0,$start1_0,$end1_0,$strand1,$identity1,$length1);
        }elsif($end1 > $end){
            if($chr1_0 eq $chr_0 && $start1_0 < $end_0 && $end1_0 > $end_0){
            my $end = $end1;
            my $end_0 = $end1_0;
            $identity = sprintf "%.4f",(($identity1*$length1+$identity*$length)/($length+$length1));
            $length = sprintf "%d",($end-$start+1);
            }else{
                if($length1*$identity1 > $length*$identity){
                    ($start,$end,$chr_0,$start_0,$end_0,$strand,$identity,$length) = ($start1,$end1,$chr1_0,$start1_0,$end1_0,$strand1,$identity1,$length1);
                }
            }
        }
    }
}
print OUT "$chr\t$start\t$end\t$chr_0\t$start_0\t$end_0\t$strand\t$identity\n";

close IN; close OUT;
}

sub selectIden{
use lib '.';
use base;

my ($mergefile,$idenfile,$iden_threshold) = @_;
my $idenfile_sort = $idenfile.".sort";
my $idenfile_merge = $idenfile.".merge";

open(IN,$mergefile)||die("seqAlign.selectIden: error with opening $mergefile\n");
open(OUT, ">$idenfile")||die("seqAlign.selectIden: error with writing to $idenfile\n");

while(<IN>){
    chomp();
    my @arr = split(/\t/,$_);
    if($arr[7] >= $iden_threshold){
	print OUT "$arr[0]\t$arr[1]\t$arr[2]\n";
    }
}
`sort -k1,1 -k2,2n $idenfile > $idenfile_sort`;
my $num=`cat $idenfile_sort |wc -l`;
if($num ne 0){
    base::merge($idenfile_sort,$idenfile_merge);
    `mv $idenfile_merge $idenfile`;
    `rm -rf $idenfile_sort`;
}
}

1;
