use strict;
use warnings;
use Getopt::Long;

my ($infile,$ems,$minDP,$minAD,$outfile,$p1,$p2,$s1,$s2,$sim);
GetOptions(
	"vcf:s"=>\$infile,
	"o:s"=>\$outfile,
	"parent1:s"=>\$p1,
	"parent2:s"=>\$p2,
	"bulk1:s"=>\$s1,
	"bulk2:s"=>\$s2,
	"sim:s"=>\$sim
);


if(!defined $infile or !defined $outfile){
my $usage=<<"USAGE";
Program: $0
Contact: BioMarble <biomarble\@163.com>
Discription:
Usage:
	Options:
	-vcf	<file>	required	main vcf file, vcf format, after filtered
	-o	<file>	required	main output file ,vcf format
	-parent1	<string>	optional	Parent2 name Mut
	-parent2	<string>	optional	Parent1 name WT
	-bulk1	<string>	optional	offspring pool 2 name Mut 
	-bulk2	<string>	optional	offspring pool 1 name WT 
	-sim	<file>	optional	simulation file 
USAGE
print $usage;
exit;
}
my %sim;
if(defined $sim){
	
	open IN,$sim or die $!;
	<IN>;
	while(<IN>){
		chomp;
		next if($_=~/^\s*$/);
		my ($dp,$str)=split /\t/,$_,2;
		$sim{$dp}=$str;
	}
	close IN;
}



open IN,$infile or die "cannot open input file \n";
open OUT,">$outfile" or die "cannot open input file\n";
if(defined $sim){
	print OUT "#CHR\tPOS\tSNPindexMut\tSNPindexWT\t△SNPindex\tL95\tH95\tL99\tH99\n";
}else{
	print OUT "#CHR\tPOS\tSNPindexMut\tSNPindexWT\t△SNPindex\n";
}
my @header;
while(<IN>){
	chomp;
	next if($_=~/^\s*$/);
	if($_=~/#/){
		if($_=~/#CHROM/){
			(undef,undef,undef,undef,undef,undef,undef,undef,undef,@header)=split /\s+/,$_;
		}
		next;
	}
	$_=~s/\|/\//g;
	my ($chr,$pos,undef,$ref,$alt,undef,undef,$str,$format,@sample)=split /\s+/,$_;
	my ($gt,$ad)=&getInfo($format,$str,@sample);
	
	my $a1=$ad->{$s1}{1-$gt->{$p1}};
	my $b1=$ad->{$s1}{$gt->{$p1}};
	my $dp1=$a1+$b1;
	my $index1=$b1/$dp1;
	my $a2=$ad->{$s2}{1-$gt->{$p1}};
	my $b2=$ad->{$s2}{$gt->{$p1}};
	my $dp2=$a2+$b2;
	my $index2=$b2/$dp2;
	my $delta=$index1-$index2;
	my $meanDP=int(($dp1+$dp2)/2);
	$meanDP=300 if($meanDP>300);
	if(defined $sim){
		print OUT "$chr\t$pos\t$index1\t$index2\t$delta\t$sim{$meanDP}\n";
	}else{
		print OUT "$chr\t$pos\t$index1\t$index2\t$delta\n";
	}
}
close IN;
close OUT;


sub getInfo{
	my ($format,$str,@sample)=@_;
	my %ad;
	if($format=~/\bAD\b/){
		for(my $i=0;$i<@sample;$i++){
			($ad{$header[$i]}{0},$ad{$header[$i]}{1})=SamAD($format,$sample[$i]);
		}
	}else{
		die "no AD format\n";
	}
	my %gt;
	for(my $i=0;$i<@sample;$i++){
		$gt{$header[$i]}=SamGT($format,$sample[$i]);
	}
	return (\%gt,\%ad);
}


sub SamGT{
	my ($format,$str)=@_;
	my @f=split /:/,$format;
	my @h=split /:/,$str;
	for(my $i=0;$i<@f;$i++){
		if($f[$i] eq "GT"){
			my $gt=(split /\//,$h[$i])[0];
			return $gt;
		}
	}
}
sub SamAD{
	my ($format,$str)=@_;
	my @f=split /:/,$format;
	my @h=split /:/,$str;
	my ($a,$b);
	for(my $i=0;$i<@f;$i++){
		if($f[$i] eq "AD"){
			
			if($h[$i] ne "."){
				($a,$b)=split /,/,$h[$i];
			}else{
				$a=0;
				$b=0;
			}
			return ($a,$b);
		}
	}
}
