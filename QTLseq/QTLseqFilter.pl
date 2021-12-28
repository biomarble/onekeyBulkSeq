
use strict;
use warnings;
use Getopt::Long;

my ($infile,$ems,$minDP1,$minDP2,$outfile,$p1,$p2,$s1,$s2,$maxDP1,$maxDP2);
GetOptions(
	"vcf:s"=>\$infile,
	"ems"=>\$ems,
	"minDP1:s"=>\$minDP1,
	"minDP2:s"=>\$minDP2,
	"maxDP1:s"=>\$maxDP1,
	"maxDP2:s"=>\$maxDP2,
	"parent1:s"=>\$p1,
	"parent2:s"=>\$p2,
	"bulk1:s"=>\$s1,
	"bulk2:s"=>\$s2,
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
	-parent1	<string>	optional	Parent1 name
	-parent2	<string>	optional	Parent2 name
	-bulk1	<string>	optional	offspring pool1  name 
	-bulk2	<string>	optional	offspring pool2  name 
	-o	<file>	required	main output file ,vcf format
	-minDP1	<float>	optional	minimun total depth for parent site, default 8;
	-maxDP1	<float>	optional	maximun total depth for parent site, default 100;
	-minDP2	<float>	optional	minimun total depth for pool site, default 20;
	-maxDP2	<float>	optional	maximun total depth for pool site, default 300;
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
	my @useheader;
	for(my $i=0;$i<@header;$i++){
		if(defined $p1 and defined $p2){
			if($header[$i] eq $s1 or $header[$i] eq $p1 or $header[$i] eq $s2 or $header[$i] eq $p2){
				push @used,$orisample[$i];
				push @useheader,$header[$i];
			}
		}else{
			if($header[$i] eq $s1 or $header[$i] eq $s2){
				push @used,$orisample[$i];
				push @useheader,$header[$i];
			}
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
	my ($gt,$ad)=&getInfo($format,$str,\@useheader,@used);
	if(defined $p1 and defined $p2){
		if($gt->{$p1} eq "0/1" or $gt->{$p2} eq "0/1"){
			$flag{'ParentHeter'}++;
			next;
		}
		if($gt->{$p1} eq $gt->{$p2} and $gt->{$p1} ne "0/1"){
			$flag{'ParentHomoEqual'}++;
			next;
		}
		if($ad->{$p1}{'ref'}+$ad->{$p1}{'alt'}>$maxDP1 or $ad->{$p2}{'ref'}+$ad->{$p2}{'alt'}>$maxDP1){
			$flag{'ParentminDP.too.high'}++;
			next;
		}
		if($ad->{$p1}{'ref'}+$ad->{$p1}{'alt'}<$minDP1 or $ad->{$p2}{'ref'}+$ad->{$p2}{'alt'}<$minDP1){
			$flag{'ParentminDP.too.low'}++;
			next;
		}
	}
	if($gt->{$s1} ne "0/1" and $gt->{$s2} ne "0/1" and $gt->{$s1}  eq $gt->{$s2} ){
		$flag{'PoolHomoEqual'}++;
		next;
	}
	if($ad->{$s1}{'ref'}+$ad->{$s1}{'alt'}>$maxDP2 or $ad->{$s2}{'ref'}+$ad->{$s2}{'alt'}>$maxDP2){
		$flag{'PoolminDP.too.high'}++;
		next;
	}
	if($ad->{$s1}{'ref'}+$ad->{$s1}{'alt'}<$minDP2 or $ad->{$s2}{'ref'}+$ad->{$s2}{'alt'}<$minDP2 ){
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
	my ($format,$str,$hea,@sample)=@_;
	my %ad;
	my @k=@{$hea};
	if($format=~/\bAD\b/){
		for(my $i=0;$i<@sample;$i++){
			($ad{$k[$i]}{'ref'},$ad{$k[$i]}{'alt'})=SamAD($format,$sample[$i]);
		}
	}else{
		die "no AD format\n";
	}
	my %gt;
	for(my $i=0;$i<@sample;$i++){
		$gt{$k[$i]}=SamGT($format,$sample[$i]);
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