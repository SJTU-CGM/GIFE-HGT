#!/usr/bin/perl
package splitDB;

sub splitDB{
use strict;
use warnings;
use base;
use Getopt::Long;

my $usage="\nUsage: GIFEHGT splitDB --id <genome_id> --taxo <genome_taxonomy> --taxoes <taxonomy_file> --type <type_file>

splitDB program is used to split the genomic database into three groups: SG, CRG and DRG.

Necessary input description:

  genome_id		<string>        The genome id of target genome.

  genome_taxonomy	<string>	The taxonomy of target genome (a comma-separated string), whose format is like [kingdom,phylum,class,order].

  taxonomy_file         <string>        The taxonomy of all genomes. each line of which is like [genome_id	genome_kingdom 	genoe_phylum	genome_class	genome_order].

  type_file		<string>	The kingdom of all genomes. Like [GCF_000001405.40_GRCh38.p14_genomic      Metazoa] for human.

Options (defaults in parentheses):

    --help                              Print this usage page.

  Output:

    --outdir		<string>	The result files will be output to this directory. (./splitDB)

  Other parameters:

    --distant		<string>	The distantly related group. kingdom or phylum can be chosen. (kingdom)

";

#read splitDB parameters
my $distant = "kingdom";
my $out_dir = "./splitDB/";
my $type_file;
my $genome_id;
my $taxanomy_target;
my $taxanomy_file;
my $help;
GetOptions(
    'id=s'=> \$genome_id,
    'taxo=s'=> \$taxanomy_target,
    'taxoes=s'=> \$taxanomy_file,
    'type=s'=> \$type_file,
    'outdir=s'  => \$out_dir,
    'help'      => \$help,
    'distant=s' => \$distant
) or die $usage;

die $usage if defined($help);
die $usage if(!$genome_id);
die $usage if(!$taxanomy_target);
die $usage if(!$taxanomy_file);
die $usage if(!$type_file);

#Check existence of output directory
if(-e $out_dir){
    print("Warning: output directory \"$out_dir\" already exists. Existing files will be overwritten.\n");
    `rm -rf $out_dir`;
}

#Adjust directory names and create output directory
unless($out_dir=~/\/$/){
    $out_dir.="/";
}
mkdir($out_dir);

my ($kingdom,$phylum,$class,$order) = split(/,/,$taxanomy_target);

#get genome id in same taxonomy with target genome
getid($type_file,$taxanomy_file,$out_dir,$genome_id,$kingdom,$phylum,$class,$order);

#get type file
gettype($type_file,$out_dir,$phylum,$class,$order);

#get self file
getself($distant,$out_dir);
}

sub getid{
my ($type_file,$taxanomy_file,$out_dir,$genome_id,$kingdom,$phylum,$class,$order)=@_;

open(TYPE,$type_file)||die("splitDB.getid: error with opening the type file\n");
open(TAXA,$taxanomy_file)||die("splitDB.getid: error with opening the taxonomy file\n");
open(KING,">".$out_dir."kingdom.id")||die("splitDB.getid: error with writing to the output directory\n");
open(EKING,">".$out_dir."e_kingdom.id")||die("splitDB.getid: error with writing to the output directory\n");
open(PHY,">".$out_dir."phylum.id")||die("splitDB.getid: error with writing to the output directory\n");
open(EPHY,">".$out_dir."e_phylum.id")||die("splitDB.getid: error with writing to the output directory\n");
open(CLA,">".$out_dir."class.id")||die("splitDB.getid: error with writing to the output directory\n");
open(ECLA,">".$out_dir."e_class.id")||die("splitDB.getid: error with writing to the output directory\n");
open(OR,">".$out_dir."order.id")||die("splitDB.getid: error with writing to the output directory\n");
open(EOR,">".$out_dir."e_order.id")||die("splitDB.getid: error with writing to the output directory\n");

@typefile = <TYPE>;
%id = ();
foreach $i (@typefile){
    chomp $i;
    @type = split(/\t/,$i);
    @name = split(/_/,$type[0]);
    $gcf = "$name[0]_$name[1]";
    $id{$gcf} = $type[0];
}

while(<TAXA>){
    chomp($_);
    my $range = $_;
    @range = split(/\t/,$range);
    if($range[0] ne $genome_id){
    if($range[1] eq $kingdom){
        print KING "$id{$range[0]}\n";
    }else{
        print EKING "$id{$range[0]}\n";
    }
    if($range[2] eq $phylum){
        print PHY "$id{$range[0]}\n";
    }else{
        print EPHY "$id{$range[0]}\n";
    }
    if($range[3] eq $class){
        print CLA "$id{$range[0]}\n";
    }else{
        print ECLA "$id{$range[0]}\n";
    }
    if($range[4] eq $order){
        print OR "$id{$range[0]}\n";
    }else{
        print EOR "$id{$range[0]}\n";
    }
    }
}

close TYPE; close TAXA; close KING; close EKING; close PHY; close EPHY; close CLA; close ECLA; close OR; close EOR;
}

sub gettype{
my ($type_file,$out_dir,$phylum,$class,$order)=@_;

open(TYPE,$type_file)||die("splitDB.gettype: error with opening the type file\n");
open(PHY,$out_dir."phylum.id")||die("splitDB.gettype: error with opening the phylum file\n");
open(CLA,$out_dir."class.id")||die("splitDB.gettype: error with opening the class file\n");
open(OR,$out_dir."order.id")||die("splitDB.gettype: error with opening the order file\n");
open(OUT,">".$out_dir."type.txt")||die("splitDB.gettype: error with writing to the output directory\n");

%type=();
while(<TYPE>){
    chomp();
    $_=~s/vertebrate_mammalian/Metazoa/;
    $_=~s/vertebrate_other/Metazoa/;
    $_=~s/invertebrate/Metazoa/;
    @arr=split(/\t/,$_);
    $type{$arr[0]}=$arr[1];
}

while(<PHY>){
    chomp();
    $type{$_}=$phylum;
}

while(<CLA>){
    chomp();
    $type{$_}=$class;
}

while(<OR>){
    chomp();
    $type{$_}=$order;
}

foreach $k (keys %type){
    print OUT "$k\t$type{$k}\n";
}

close TYPE; close PHY; close CLA; close OR; close OUT;
}

sub getself{
use lib '.';
use base;
my ($distant,$out_dir)=@_;
my $kingdomid=$out_dir."kingdom.id";
my $phylumid=$out_dir."phylum.id";
my $classid=$out_dir."class.id";
my $orderid=$out_dir."order.id";

if($distant eq "kingdom"){
    my $kingdomSpeciesid=$out_dir."kingdom-species.id";
    my $kingdomOrderid=$out_dir."kingdom-order.id";
    my $kingdomClassid=$out_dir."kingdom-class.id";
    my $kingdomPhylumid=$out_dir."kingdom-phylum.id";
    my $comcp = "cp $kingdomid $kingdomSpeciesid";
    system($comcp);
    base::notin($orderid,$kingdomid,$kingdomOrderid);
    base::notin($classid,$kingdomid,$kingdomClassid);
    base::notin($phylumid,$kingdomid,$kingdomPhylumid);
}elsif($distant eq "phylum"){
    my $phylumSpeciesid=$out_dir."phylum-species.id";
    my $phylumOrderid=$out_dir."phylum-order.id";
    my $phylumClassid=$out_dir."phylum-class.id";
    my $comcp = "cp $phylumid $phylumSpeciesid";
    system($comcp);
    base::notin($orderid,$phylumid,$phylumOrderid);
    base::notin($classid,$phylumid,$phylumClassid);
}elsif($distant eq "class"){
    my $classSpeciesid=$out_dir."class-species.id";
    my $classOrderid=$out_dir."class-order.id";
    my $comcp = "cp $classid $classSpeciesid";
    system($comcp);
    base::notin($orderid,$classid,$classOrderid);
}else{
    print "Wrong parameter for --distant. Only kingdom and phylum can be chosen.";
}
}

1;
