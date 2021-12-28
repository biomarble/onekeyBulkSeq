use strict;
use warnings;
use Getopt::Long;

my ($infile,$ems,$minDP,$minAD,$outfile,$power,$s1,$s2);
GetOptions(
	"vcf:s"=>\$infile,
	"o:s"=>\$outfile,
	"power:s"=>\$power,
	"bulk1:s"=>\$s1,
	"bulk2:s"=>\$s2,
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
	-power	<string>	optional	power default 4
	-bulk1	<string>	optional	offspring pool 2 name Mut 
	-bulk2	<string>	optional	offspring pool 1 name WT 
USAGE
print $usage;
exit;
}

$power=4 if(!defined $power);

open IN,$infile or die "cannot open input file \n";
open OUT,">$outfile" or die "cannot open input file\n";
print OUT "#CHR\tPOS\tED\tED$power\n";
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
	
	my $a1=$ad->{$s1}{0};
	my $b1=$ad->{$s1}{1};
	my $to1=$a1+$b1;

	my $a2=$ad->{$s2}{0};
	my $b2=$ad->{$s2}{1};
	my $to2=$a2+$b2;

	my $ed=sqrt(($a1/$to1-$a2/$to2)**2+($b1/$to1-$b2/$to2)**2);
	my $ed2=$ed**$power;

	print OUT "$chr\t$pos\t$ed\t$ed2\n";
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
			
			if($h[$i] ne "." and $h[$i] ne './.' ){
				($a,$b)=split /,/,$h[$i];
			}else{
				$a=0;
				$b=0;
			}
			return ($a,$b);
		}
	}
}
