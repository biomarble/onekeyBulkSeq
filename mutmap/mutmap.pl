#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);

my $version="1.0.0";
# ==============================================================
# Get Options
# ==============================================================
my ($vcf,$ref,$ems,$minDP1,$minDP2,$od,$win,$step,$mode,$key,$thr,$minCount,$maxDP1,$maxDP2,$p,$s,);
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
	"key:s"=>\$key,
	"minDP1:s"=>\$minDP1,
	"minDP2:s"=>\$minDP2,
	"maxDP1:s"=>\$maxDP1,
	"maxDP2:s"=>\$maxDP2,
	"minCount:s"=>\$minCount,
	"parent:s"=>\$p,
	"bulk:s"=>\$s,
	"ems"=>\$ems,
	"win:s"=>\$win,
	"step:s"=>\$step,
	"thr:s"=>\$thr,
	"od:s"=>\$od,
	) or &USAGE;
&USAGE unless ($vcf and $p and $s);

$thr="NULL" if(!defined $thr);
$minDP1=8 if(!defined $minDP1);
$minDP2=20 if(!defined $minDP2);
$maxDP1=100 if(!defined $maxDP1);
$maxDP2=300 if(!defined $maxDP2);
$minCount=5 if(!defined $minCount);
$win=2000000 if(!defined $win);
$step=100000 if(!defined $step);
$key="mutmap" if(!defined $key);
$od="MutmapOut" if(!defined $od);
$od=abs_path($od);
$vcf=abs_path($vcf);
mkdir $od if(!-e $od);
my $cmd="vcftools --vcf $vcf  --recode --indv $p --indv $s --stdout >$od/0.$key.use.vcf 2>$od/0.$key.use.vcf.log && ";
$cmd.="perl $Bin/MutmapFilter.pl -vcf $od/0.$key.use.vcf -parent $p -bulk $s -maxDP1 $maxDP1 -maxDP2 $maxDP2 -minDP1 $minDP1 -minDP2 $minDP2 -o $od/0.$key.filtered.vcf";
$cmd.=" -ems" if(defined $ems);
runCMD($cmd,'SNP filteration');

$cmd="perl $Bin/MutmapIndex.pl  -vcf  $od/0.$key.filtered.vcf -o $od/1.$key.raw.index.txt -p $p -s $s";
runCMD($cmd,'SNPindex calculation');

$cmd="Rscript $Bin/MutmapSliding.r $od/1.$key.raw.index.txt $win $step $minCount $od/2.$key.windowed.index.txt";
runCMD($cmd,'SNPindex sliding');

$cmd="Rscript $Bin/MutmapPlot.r  $od/2.$key.windowed.index.txt $od/1.$key.raw.index.txt  $thr  \"$key SNP-index Plot\"  $od/3.$key.manhattan";
runCMD($cmd,'SNPindex plot');

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
        -parent          <string>    required, 野生亲本样品名（对应VCF文件） 
        -bulk          <string>    required, 突变混池样品名（对应VCF文件） 
        -key        <string>    optional, 项目关键词，default:key
        -od         <string>    optional, 输出目录，default:MutmapOut

        -minDP1     <int>       optional, 亲本最低深度，default:8
        -minDP2     <int>       optional, 子代混池最低深度，default:20
        -maxDP1     <int>       optional, 亲本最高深度, default:100
        -maxDP2     <int>       optional, 子代混池最高深度，default 300
        -minCount   <int>       optional, 拟合窗口最少SNP数目，default 5

        -win        <int>       optional, 拟合窗口大小, default 2000000
        -step       <int>       optional, 拟合步长，default 100000
        -ems        <option>    optional，是否EMS诱变

        -thr        <float>     optional，设置阈值，用于绘图，默认不绘制
        -help       help
USAGE
                    print $usage;
                    exit;
}
