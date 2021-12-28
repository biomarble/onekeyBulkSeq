use strict;
use warnings;
use Getopt::Long;

my ($infile,$ems,$minDP,$minAD,$outfile,$p1,$p,$s,$s2);
GetOptions(
	"vcf:s"=>\$infile,
	"o:s"=>\$outfile,
	"p:s"=>\$p,
	"s:s"=>\$s,
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
	-p2	<string>	optional	Parent2 name WT
	-s1	<string>	optional	offspring pool 1 name Mut 
USAGE
print $usage;
exit;
}

open IN,$infile or die "cannot open input file \n";
open OUT,">$outfile" or die "cannot open input file\n";
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
	
	my $a=$ad->{$s}{$gt->{$p}};
	my $b=$ad->{$s}{1-$gt->{$p}};
	my $index=$b/($a+$b);
	print OUT "$chr\t$pos\t$index\n";
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