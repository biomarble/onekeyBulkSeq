#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);

my $version="1.0.0";
my ($vcf,$ref,$ems,$minDP1,$minDP2,$od,$win,$step,$mode,$key,$thr,$minCount,$maxDP1,$maxDP2,$p1,$s1,$p2,$s2,$method,$power,$repnum,$pop,$nPool,$thr1,$thr2);
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
	"key:s"=>\$key,
	"minDP1:s"=>\$minDP1,
	"minDP2:s"=>\$minDP2,
	"maxDP1:s"=>\$maxDP1,
	"maxDP2:s"=>\$maxDP2,
	"minCount:s"=>\$minCount,
	"parent1:s"=>\$p1,
	"parent2:s"=>\$p2,
	"bulk1:s"=>\$s1,
	"bulk2:s"=>\$s2,
	"win:s"=>\$win,
	"step:s"=>\$step,
	"thr:s"=>\$thr,
	"method:s"=>\$method,
	"plotthrED:s"=>\$thr1,
	"plotthrIndex:s"=>\$thr2,
	"power:s"=>\$power,
	"od:s"=>\$od,
	) or &USAGE;
&USAGE unless ($vcf  and $s1 and $s2);

$method="both" if(!defined $method);
$minDP1=8 if(!defined $minDP1);
$minDP2=20 if(!defined $minDP2);
$maxDP1=100 if(!defined $maxDP1);
$maxDP2=300 if(!defined $maxDP2);
$minCount=5 if(!defined $minCount);
$win=2000000 if(!defined $win);
$step=100000 if(!defined $step);
$power=4 if(!defined $power);
$key="BSA" if(!defined $key);
$od="BSAOut" if(!defined $od);
$od=abs_path($od);
$vcf=abs_path($vcf);
$repnum=1000 if(!defined $repnum);
$nPool=40 if(!defined $nPool);
$pop='F2' if(!defined $pop);
$thr1="NULL" if(!defined $thr1);
$thr2="NULL" if(!defined $thr2);

mkdir $od if(!-e $od);
my $cmd;

if(defined $p1 and defined $p2){
  $cmd="vcftools --vcf $vcf  --recode --indv $p1 --indv $p2 --indv $s1 --indv $s2 --stdout >$od/0.$key.use.vcf 2>$od/0.$key.use.vcf.log && ";
  $cmd.="perl $Bin/QTLseqFilter.pl -vcf $od/0.$key.use.vcf -parent1 $p1 -parent2 $p2 -bulk1 $s1 -bulk2 $s2 -maxDP1 $maxDP1 -maxDP2 $maxDP2 -minDP1 $minDP1 -minDP2 $minDP2 -o $od/0.$key.filtered.vcf";
}else{
  $cmd="vcftools --vcf $vcf  --recode --indv $s1 --indv $s2 --stdout >$od/0.$key.use.vcf 2>$od/0.$key.use.vcf.log && ";
  $cmd.="perl $Bin/QTLseqFilter.pl -vcf $od/0.$key.use.vcf -bulk1 $s1 -bulk2 $s2 -maxDP2 $maxDP2 -minDP2 $minDP2 -o $od/0.$key.filtered.vcf";
}
$cmd.=" -ems" if(defined $ems);

runCMD($cmd,'SNP filteration');

if($method eq "index" or $method eq "both"){
if(defined $p1 and defined $p2){
  if($thr2 eq "NULL"){
    $cmd="Rscript $Bin/QTLseqIndexSim.r $nPool $repnum 0 $pop $od/0.$key.simulation.index.txt";
    runCMD($cmd,'SNPindex simulation');
    $cmd="perl $Bin/QTLseqIndex.pl  -vcf  $od/0.$key.filtered.vcf -o $od/1.$key.raw.index.txt -parent1 $p1 -parent2 $p2 -bulk1 $s1 -bulk2 $s2   -sim $od/0.$key.simulation.index.txt";
    runCMD($cmd,'SNPindex calculation');
    $cmd="Rscript $Bin/QTLseqSlidingIndex.r $od/1.$key.raw.index.txt $win $step $minCount $od/2.$key.windowed.index.txt  1";
    runCMD($cmd,'SNPindex sliding');
  }else{
    $cmd="perl $Bin/QTLseqIndex.pl  -vcf  $od/0.$key.filtered.vcf -o $od/1.$key.raw.index.txt -parent1 $p1 -parent2 $p2 -bulk1 $s1 -bulk2 $s2";
    runCMD($cmd,'SNPindex calculation');
    $cmd="Rscript $Bin/QTLseqSlidingIndex.r $od/1.$key.raw.index.txt $win $step $minCount $od/2.$key.windowed.index.txt  0";
    runCMD($cmd,'SNPindex sliding');
  }
  $cmd="Rscript $Bin/QTLseqPlotIndex.r  $od/2.$key.windowed.index.txt $od/1.$key.raw.index.txt  $thr2  \"$key SNP-index \"  $od/3.$key.index.manhattan "; 
  runCMD($cmd,'SNPindex plot');
}else{
  warn("Cannot run SNP-index without parents!\n");
}
}

if($method eq "ed" or $method eq "both"){
  $cmd="perl $Bin/QTLseqED.pl  -vcf  $od/0.$key.filtered.vcf -o $od/1.$key.raw.ED.txt -bulk1 $s1 -bulk2 $s2 -power $power";
  runCMD($cmd,'ED calculation');
  $cmd="Rscript $Bin/QTLseqSlidingED.r $od/1.$key.raw.ED.txt $win $step $minCount $od/2.$key.windowed.ED.txt";
  runCMD($cmd,'ED sliding');
  $cmd="Rscript $Bin/QTLseqPlotED.r  $od/2.$key.windowed.ED.txt $od/1.$key.raw.ED.txt  $thr1  \"$key Euclidean Distance plot (Power=$power)\"  $od/3.$key.ED.manhattan";
  runCMD($cmd,'ED plot');
}


sub runCMD{
	my ($cmd,$log)=@_;
	print "$log is running...\n";
	my $flag=system $cmd;
	if($flag!=0){
		die "FATAL RUN ERROR:\n$cmd\n ";
	}
	print "$log is done!\n\n";

}

sub USAGE {
	my $usage=<<"USAGE";
	Program: $0
	Version: $version
	Contact: BioMarble <biomarble\@163.com>
	Discription:
	Usage:
        Options: 
        -vcf        <file>      required, SNP文件,vcf格式
        -parent1    <string>    required, 突变性状亲本样品名（对应VCF文件）
        -parent2    <string>    required, 野生性状亲本样品名（对应VCF文件）
        
        -bulk1      <string>    required, 突变性状混池样品名（对应VCF文件）
        -bulk2      <string>    required, 野生性状混池样品名（对应VCF文件）
        
        -key        <string>    optional, 项目关键词，default:BSA
        -od         <string>    optional, 输出目录，default:BSAOut

        -minDP1     <int>       optional, 亲本最低深度，default:8
        -minDP2     <int>       optional, 子代混池最低深度，default:20
        -maxDP1     <int>       optional, 亲本最高深度, default:100
        -maxDP2     <int>       optional, 子代混池最高深度，default 300
        -minCount   <int>       optional, 拟合窗口最少SNP数目，default 5

        -win        <int>       optional, 拟合窗口大小, default 2000000
        -step       <int>       optional, 拟合步长，default 100000
        -power      <int>       optional, ED power，default 4

        -plotthrED           <float>     optional，用于ED绘图，默认为 median+3*SD
        -plotthrIndex        <float>     optional，用于SNP-index绘图，默认为计算机模拟实验99%与95%置信区间
        -help       help
USAGE
                    print $usage;
                    exit;
}
