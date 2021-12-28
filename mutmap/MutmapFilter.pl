#perl $0 -ems -dp 4 -ad 2 
#input vcf  , output vcf
use strict;
use warnings;
use Getopt::Long;

my ($infile,$ems,$minDP1,$minDP2,$outfile,$p,$p2,$s,$s2,$maxDP1,$maxDP2);
GetOptions(
	"vcf:s"=>\$infile,
	"ems"=>\$ems,
	"minDP1:s"=>\$minDP1,
	"minDP2:s"=>\$minDP2,
	"maxDP1:s"=>\$maxDP1,
	"maxDP2:s"=>\$maxDP2,
	"parent:s"=>\$p,
	"bulk:s"=>\$s,
	"o:s"=>\$outfile,
);
if(!defined $infile or !defined $outfile){
my $usage=<<"USAGE";
Program: $0
Contact: BioMarble <biomarble\@163.com>
Discription:
Usage:
	Options:
	-vcf	<file>	required	main vcf file, vcf format, 
	-o	<file>	required	main output file ,vcf format
	-ems	<opt>	optional	whether to filter SNPs by EMS criteria( C/G to T/A )
	-minDP1	<float>	optional	minimun total depth for parent site, default 8;
	-minDP2	<float>	optional	minimun total depth for pool site, default 20;
	-maxDP1	<float>	optional	maximun total depth for pool site, default 100;
	-maxDP2	<float>	optional	maximun total depth for pool site, default 300;
	-p	<string>	optional	Parent name
	-s	<string>	optional	offspring pool  name 
USAGE
print $usage;
exit;
}

open IN,$infile or die "cannot open input file: $infile \n";
open OUT,">$outfile" or die "cannot write output file: $outfile\n";

my @header;
my %flag;

$minDP1=8 if(!defined $minDP1);
$minDP2=20 if(!defined $minDP2);
$maxDP1=100 if(!defined $maxDP1);
$maxDP2=300 if(!defined $maxDP2);
my $rem;

while(<IN>){
	chomp;
	next if($_=~/^\s*$/);
	$_=~s/\|/\//g;
	if($_=~/#/){
		if($_=~/#CHROM/){
			(undef,undef,undef,undef,undef,undef,undef,undef,undef,@header)=split /\s+/,$_;
		}
		print OUT $_."\n";
		next;
	}
	my ($chr,$pos,undef,$ref,$alt,undef,undef,$str,$format,@orisample)=split /\s+/,$_;
	if($alt=~/,/){
		$flag{"multAlt"}++;
		next;
	}
	my @used;
	for(my $i=0;$i<@header;$i++){
		if($header[$i] eq $s or $header[$i] eq $p){
			push @used,$orisample[$i];
		}
	}
	
	if(join(",",@used) =~/\.\/\./ ){
		$flag{"missingGT"}++;
		next;
	}
	
	if(defined $ems){
		if($ref eq "G" and $alt eq "A" or $ref eq "C" and $alt eq "T" or $ref eq "A" and $alt eq "G" or $ref eq "T" and $alt eq "C"){
		}else{
			$flag{"EMS"}++;
			next;
		}
	}
	my ($gt,$ad)=&getInfo($format,$str,@orisample);
	
	if($gt->{$p} eq "0/1"){
		$flag{'ParentHeter'}++;
		next;
	}
	if($gt->{$p} eq $gt->{$s} and $gt->{$p} ne "0/1"){
		$flag{'ParentPoolHomoEqual'}++;
		next;
	}
	if($ad->{$p}{'ref'}+$ad->{$p}{'alt'}>$maxDP1 ){
		$flag{'ParentminDP.too.high'}++;
		next;
	}
	if($ad->{$p}{'ref'}+$ad->{$p}{'alt'}<$minDP1){
		$flag{'ParentminDP.too.low'}++;
		next;
	}
	
	if($ad->{$s}{'ref'}+$ad->{$s}{'alt'}>$maxDP2 ){
		$flag{'PoolminDP.too.high'}++;
		next;
	}
	if($ad->{$s}{'ref'}+$ad->{$s}{'alt'}<$minDP2 ){
		$flag{'PoolminDP.too.low'}++;
		next;
	}
	$rem++;
	print OUT "$_\n";
}
close OUT;

open LOG,">$outfile.filter.log" or die "cannot write log file $outfile.filter.log\n";
foreach my $key (sort keys %flag){
	print LOG "$key:\t$flag{$key}\n";
}
print LOG "finalSNP:\t$rem\n";
close LOG;





sub getInfo{
	my ($format,$str,@sample)=@_;
	my %ad;
	if($format=~/\bAD\b/){
		for(my $i=0;$i<@sample;$i++){
			($ad{$header[$i]}{'ref'},$ad{$header[$i]}{'alt'})=SamAD($format,$sample[$i]);
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
			return $h[$i];
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