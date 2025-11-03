#!/usr/bin/perl
package base;

#Get fasta format file from bed format file
sub getSeq{
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

#Get lines in file2 but not in file1
sub notin{
my ($id1,$id2,$out)= @_;

my %hash = ();

open(ID1,$id1)||die("notin: error with opening $id1\n");
open(ID2,$id2)||die("notin: error with opening $id2\n");
open(OUT,">$out")||die("notin: error with writing $out\n");
while(<ID1>){
    chomp();
    $hash{$_} = 1;
}

while(<ID2>){
    chomp();
    if(not exists($hash{$_})){
        print OUT "$_\n";
    }
}
close ID1;close ID2; close OUT;
}

#Get lines in file2 and file1
sub in{
my ($id1,$id2,$out)= @_;

my %hash = ();

open(ID1,$id1)||die("notin: error with opening $id1\n");
open(ID2,$id2)||die("notin: error with opening $id2\n");
open(OUT,">$out")||die("notin: error with writing $out\n");
while(<ID1>){
    chomp();
    $hash{$_} = 1;
}

while(<ID2>){
    chomp();
    if(exists($hash{$_})){
        print OUT "$_\n";
    }
}
close ID1;close ID2; close OUT;
}

#Get info file2 of some id in file1
sub id2info{
my ($id1,$info2,$out)= @_;

my %hash = ();

open(ID1,$id1)||die("notin: error with opening $id1\n");
open(INFO2,$info2)||die("notin: error with opening $id2\n");
open(OUT,">$out")||die("notin: error with writing $out\n");
while(<INFO2>){
    chomp();
    my @arr = split(/\t/,$_);
    $hash{$arr[0]} = $_;
}

while(<ID1>){
    chomp();
    if(exists $hash{$_}){
        print OUT "$hash{$_}\n";
    }
}
close ID1;close INFO2; close OUT;
}


#Sort file in bed or id format
sub sort{
my ($in,$out) = @_;

# 检查输入文件是否为空（如果不存在也视为空）
if (! -e $in || -z $in) {
    open(OUT, ">", $out) or die "sort: cannot write to $out: $!\n";
    close OUT;  # 创建一个空文件
    return;     # 直接返回，不执行后续操作
}

my $in_re = $in.".replace";
open(INR,">$in_re")||die("sort: error with writing $in_re\n");
open(IN,$in)||die("sort: error with opening $in\n");
while(<IN>){
    chomp($_);
    my @arr = split(/\t|-/,$_);
    print INR "$arr[0]\t$arr[1]\t$arr[2]\n";
}
close INR;

my $comsort = "sort -k1,1 -k2,2n $in_re > $out";
system($comsort);
my $comrm = "rm -rf $in_re";
system($comrm);
}

#Merge file in bed format
sub merge{
my ($in,$out) = @_;
open(IN,$in)||die("merge: error with opening $in\n");
open(OUT,">$out")||die("merge: error with writing $out\n");

my ($chr,$start,$end) = ("",-1,-1);
my $has_data = 0;  # 标记是否有有效数据

while(<IN>){
    chomp($_);
    $has_data = 1;  # 至少有一行数据
    my ($chr1,$start1,$end1) = split(/\s+/,$_);
    if($chr eq ""){
        ($chr,$start,$end) = ($chr1,$start1,$end1);
    }
    elsif($chr1 ne $chr){
        print OUT "$chr\t$start\t$end\n";
        ($chr,$start,$end) = ($chr1,$start1,$end1);
    }
    else{
        if($start1 > $end +1){
            print OUT "$chr\t$start\t$end\n";
            ($start,$end) = ($start1,$end1);
        }
        elsif($end1 > $end){
            $end = $end1;
        }
    }
}
print OUT "$chr\t$start\t$end\n" if $has_data;
close IN;close OUT;
}

#Get fasta information format file from bed format file
sub getHGTseqInfo{
my ($fa,$bed,$out)= @_;
open(FA,$fa)||die("bed2fasta: error with opening $fa\n");
open(BED,$bed)||die("bed2fasta: error with opening $bed\n");
open(OUT,">$out")||die("bed2fasta: error with writing to $out\n");

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
    my ($chr,$start,$end,$info) = ($arr[0],$arr[1],$arr[2],$arr[3]);
    my $HGT = uc(substr($hash{$chr},$start-1,$end-$start+1));
    print OUT ">$chr-$start-$end\t$info\n$HGT\n";
}
close FA;close BED;close OUT;
}

#Get bed format file from fasta format file
sub fastatobed{
my ($in,$out) = @_;
open(IN,$in)||die("fasta2bed: error with opening $in\n");
open(OUT,">$out")||die("fasta2bed: error with writing to $out\n");

while(<IN>){
    chomp($_);
    if($_ =~ />(.*)/){
	my @arr = split(/-|:| /,$1);
        print OUT "$arr[0]\t$arr[1]\t$arr[2]\n";
    }
}
close IN; close OUT;
}

#Convert bed format to id format
sub bed2id{
my ($in,$out) = @_;
open(IN,$in)||die("bed2id: error with opening $in\n");
open(OUT,">$out")||die("bed2id: error with writing to $out\n");

while(<IN>){
    chomp($_);
    my @arr = split(/\t/,$_);
    print OUT "$arr[0]-$arr[1]-$arr[2]\n";
}
close IN; close OUT;
}

#Convert id format to bed format
sub id2bed{
my ($in,$out) = @_;
open(IN,$in)||die("id2bed: error with opening $in\n");
open(OUT,">$out")||die("id2bed: error with writing to $out\n");

while(<IN>){
    chomp($_);
    my @arr = split(/-|\t/,$_);
    print OUT "$arr[0]\t$arr[1]\t$arr[2]\n";
}
close IN; close OUT;
}

#Get first two columns of file
sub get2col{
my ($in,$out) = @_;
open(IN,$in)||die("get2col: error with opening $in\n");
open(OUT,">$out")||die("get2col: error with writing to $out\n");

while(<IN>){
    chomp($_);
    my @arr = split(/\t| /,$_);
    print OUT "$arr[0]\t$arr[1]\n";
}
close IN; close OUT;
}

#Filter simple or low
sub filterSimple{
my ($in,$out) = @_;
open(IN,$in)||die("filterSimple: error with opening $in\n");
open(OUT,">$out")||die("filterSimple: error with writing $out\n");
while(<IN>){
    chomp();
    @arr = split(/\t/,$_);
    if(($_ =~ /Simple/ || $_ =~ /Low/) && $arr[0] ne ""){
        print OUT "$arr[0]\t$arr[1]\t$arr[2]\n";
    }
}
close IN;close OUT;
}

#compare two bed files, used to get the overlapped regions
sub biodiff{
my ($bed1,$bed2,$out)= @_;
my $bed1sort = $bed1.".sort";
my $bed2sort = $bed2.".sort";
my $comsort1 = "sort -k1,1 -k2,2n $bed1 > $bed1sort";
my $commv1 = "mv $bed1sort $bed1";
system($comsort1);
system($commv1);
my $comsort2 = "sort -k1,1 -k2,2n $bed2 > $bed2sort";
my $commv2 = "mv $bed2sort $bed2";
system($comsort2);
system($commv2);

open(BED1,$bed1)||die("biodiff: error with opening $bed1\n");
open(BED2,$bed2)||die("biodiff: error with opening $bed2\n");
open(OUT,">$out")||die("biodiff: error with writing to $out\n");
my @data = ();
my $num = 0;
my $size = 0;
while(<BED1>){
    chomp();
    my @arr = split(/\s+/,$_);
    $size = @arr;
    for(my $j=0;$j<@arr;$j++){
        $data[$num][$j] = $arr[$j];
    }
    $num += 1;
    next;
}
my $index = 0;
while(<BED2>){
    chomp();
    my @arr = split(/\s+/,$_);
    my $index_update = 0; ## whether the index is updated after this record;
    my ($chr,$start,$end) = ($arr[0],$arr[1],$arr[2]);
    for(my $i=$index;$i<$num;$i++){
        if($chr lt $data[$i][0] || ($chr eq $data[$i][0] && $end <= $data[$i][1])){
            if($index_update == 0){
                $index = $i;
            }
            last;
        }
        elsif($chr eq $data[$i][0] && (!($start >= $data[$i][2]))){ #overlap
            print OUT "$_ ; ";
            for(my $j=0;$j<$size;$j++){
                if($j==$size-1){
                    print OUT "$data[$i][$j]\n";
                }
                else{
                    print OUT "$data[$i][$j]\t";
                }
            }
            if($index_update == 0){
                $index = $i; #update the index position
                $index_update = 1;
            }
        }
    }
}
}
1;
