#/usr/bin/perl
use List::Util qw(min max);
use Statistics::Basic qw(:all);
##aligned to wild type cct
## use R,library(limma),library(vsn)
## remove the first and last 2 time points
## use glm to find the predicted value for real last time point
## use quantile(na.omit(data))[4] (75%) as 1. >1=expression
## use mean canberra distance to describe the similarity

if(@ARGV<3){die "Usage: perl $0 <input CD file folders> <input records> <output prefix>\n";}
$cellpath="/media/sf_G_DRIVE/cell-lineage/branch_cell/all_cell";
$wtccpath="/media/sf_G_DRIVE/cell-lineage/branch_cell/wt_ave_cct_v3.txt";
$treepath="/media/sf_G_DRIVE/cell-lineage/branch_cell/parent_daughters.txt";
$wtexppath="/media/sf_G_DRIVE/cell-lineage/branch_cell/wt_ave_lmexp_v3.txt";

open(LOG,">log.$$") or die "can not open log.$$,$!\n";

# define CD file folders
my $dir=$ARGV[0];
opendir DIR, $dir or die "cannot open dir $dir: $!";
my @file= readdir DIR;
closedir DIR;

# define cells and lineages
our $title="cell\tWT";
our @abaltree=();
our @abartree=();
our @abpltree=();
our @abprtree=();
our @eatree=();
our @eptree=();
our @msatree=();
our @msptree=();
our @catree=();
our @cptree=();
our @dtree=();
our @all=();
our %treearray=();
open(IN,$cellpath) or die "cannot open $cellpath: $!";
while(<IN>){
	chomp;
	@t=split(/\t/);
	if($t[0]=~/^ABal/){
		push(@abaltree,$t[0]);
		$treearray{'ABal'}.=$t[0].";";
		}
	if($t[0]=~/^ABar/){
		push(@abartree,$t[0]);
		$treearray{'ABar'}.=$t[0].";";
		}
	if($t[0]=~/^ABpl/){
		push(@abpltree,$t[0]);
		$treearray{'ABpl'}.=$t[0].";";
		}
	if($t[0]=~/^ABpr/){
		push(@abprtree,$t[0]);
		$treearray{'ABpr'}.=$t[0].";";
		}
	if($t[0]=~/^Ea/){
		push(@eatree,$t[0]);
		$treearray{'Ea'}.=$t[0].";";
		}
	if($t[0]=~/^Ep/){
		push(@eptree,$t[0]);
		$treearray{'Ep'}.=$t[0].";";
		}
	if($t[0]=~/^MSa/){
		push(@msatree,$t[0]);
		$treearray{'MSa'}.=$t[0].";";
		}
	if($t[0]=~/^MSp/){
		push(@msptree,$t[0]);
		$treearray{'MSp'}.=$t[0].";";
		}
	if($t[0]=~/^Ca/){
		push(@catree,$t[0]);
		$treearray{'Ca'}.=$t[0].";";
		}
	if($t[0]=~/^Cp/){
		push(@cptree,$t[0]);
		$treearray{'Cp'}.=$t[0].";";
		}
	if($t[0]=~/^D/){
		push(@dtree,$t[0]);
		$treearray{'D'}.=$t[0].";";
		}
	}
close(IN);

@ctree=(@catree,@cptree,'C');
@mstree=(@msatree,@msptree,'MS');
@etree=(@eatree,@eptree,'E');
@p3tree=(@dtree,'P4','P3');
@all=(@abaltree,@abartree,@abpltree,@abprtree,@ctree,@etree,@mstree,@p3tree);
@e2=('Ea','Ep');
@e4=('Eal','Ear','Epl','Epr');
@e8=('Eala','Ealp','Eara','Earp','Epla','Eplp','Epra','Eprp');
foreach $v(@e2){$e2flag{$v}=1;}
foreach $v(@e4){$e4flag{$v}=1;}
foreach $v(@e8){$e8flag{$v}=1;}
$treearray{'C'}=$treearray{'Ca'}.$treearray{'Cp'}."C";
$treearray{'E'}=$treearray{'Ea'}.$treearray{'Ep'}."E";
$treearray{'MS'}=$treearray{'MSa'}.$treearray{'MSp'}."MS";
$treearray{'P3'}=$treearray{'D'}."P3;P4";

# read in wild type cell cycle
our %wtcc=();
our %birth=();
our @e2cc=();	
our @e4cc=();	
our @e8cc=();
open(IN,$wtccpath) or die "cannot open ave cct file $wtccpath: $!";
while(<IN>){
	chomp;
	@t=split(/\t/);
	$wtcc{$t[0]}=$t[3]; ## it's in the 4th column.
	$birth{$t[0]}=$t[1];
	if(exists $e2flag{$t[0]}){push(@e2cc,$t[3]);}
	if(exists $e4flag{$t[0]}){push(@e4cc,$t[3]);}
	if(exists $e8flag{$t[0]}){push(@e8cc,$t[3]);}
	}
close(IN);
$e2ave=mean(@e2cc);
$e4ave=mean(@e4cc);
$e8ave=mean(@e8cc);

# read in cell tree relationships
%mother=();
open(IN,$treepath) or die "cannot open file $treepath: $!";
while(<IN>){
	chomp;
	@t=split(/\t/);
	$mother{$t[1]}=$t[0];
	$mother{$t[2]}=$t[0];
	}
close(IN);

# read in records
our %tpend=();
our %endcell=();
our %ko=();
our %marker=();
open(IN,$ARGV[1]) or die "cannot open file $ARGV[1]: $!";
print STDERR "reading $ARGV[1]\n";
while(<IN>){
	chomp;
	@t=split(/\t/);
	$id="CD".$t[0].".csv";
	$tpend{$id}=$t[12];
	$endcell{$id}=$t[11];
	$ko{$id}=$t[3];
	$marker{$id}=$t[2];
	}
close(IN);

# read in wild type expression
our %wexppha4=();
our %wexpnhr25=();
open(IN,$wtexppath) or die "cannot open ave exp file $wtexppath: $!";
while(<IN>){
	chomp;
	@t=split(/\t/);
	$wexppha4{$t[0]}=$t[3];
	$wexpnhr25{$t[0]}=$t[4];
	}
close(IN);

# define the CD files
my $dir=$ARGV[0];
opendir DIR, $dir or die "cannot open dir $dir: $!";
my @file= readdir DIR;
closedir DIR;

$cccor=$ARGV[2].".cc.cor.2WT";
$ccrep=$ARGV[2].".cc.cor.rep";
$expcor=$ARGV[2].".exp.cor.2WT";
$exprrep=$ARGV[2].".exp.cor.rep";
$mat1=$ARGV[2].".pha4.rawexp.txt";
$mat2=$ARGV[2].".pha4.vsnexp.txt";
$mat3=$ARGV[2].".nhr25.rawexp.txt";
$mat4=$ARGV[2].".nhr25.vsnexp.txt";
$mat5=$ARGV[2].".rawcc.txt";
$bigv=$ARGV[2].".bigVar.list.txt";
open(FILE,">$cccor") or die "can not open $cccor,$!\n";
print FILE "file\tmanually end time point\tcells edited\tKO gene\tABal\tABar\tABpl\tABpr\tC\tCa\tCp\tE\tEa\tEp\tMS\tMSa\tMSp\tD\tP3\twhole\tVariedSmallTree\tVariedBigTree\tStableSmallTree\tStableBigTree\t";
print FILE "E~MS\tC~MS\tC~AB\tAB\t";
print FILE "E2-WT\tE4-WT\tE8-WT\t";
print FILE "Cap-Caa\tCaap-Caaa\tP4\n";
close(FILE);
open(EXP,">$expcor") or die "can not open $expcor,$!\n";
print EXP "file\tmanually end time point\tcells edited\tKO gene\tABal\tABar\tABpl\tABpr\tC\tCa\tCp\tE\tEa\tEp\tMS\tMSa\tMSp\tD\tP3\twhole\tVariedSmallTree\tVariedBigTree\tStableSmallTree\tStableBigTree";
close(EXP);

open(PHA4MA,">$mat1") or die "can not open $mat1,$!\n";
open(PHA4MB,">$mat2") or die "can not open $mat2,$!\n";
open(NHR25MA,">$mat3") or die "can not open $mat3,$!\n";
open(NHR25MB,">$mat4") or die "can not open $mat4,$!\n";
open(CCM,">$mat5") or die "can not open $mat5,$!\n";

open(BV,">$bigv") or die "can not open $bigv,$!\n";
print BV "RNAi\tExpression raised\tExpression lost\n";

our $cctitle="cell\tWTcc";
our $pha4title="cell\tWTexp";
our $nhr25title="cell\tWTexp";
our %pha4exp=();
our %nhr25exp=();
our %vexp=();
our %vcc=();
$dir=~s/\/$//;
$fc=0;
our %exprep=();
our %exprepname=();
our $cor=1.5;
our $nhr25marker=0;
our $pha4marker=0;
our $expname="";
foreach $f(@file){
	if($f ne '.' && $f ne ".."){
		$cctitle.="\t".$f;
		# $expname=substr($f,0,-4);
		$expname=substr($f,0,-7);
		$expname=~s/^CD//;
		if(not exists $exprep{$expname}){$exprep{$expname}=$fc;}
		else{$exprep{$expname}.=";".$fc;}
		$exprepname{$fc}=$f;
		$fc++;
		%wtexp=();
		$ifexp=0;
		if($marker{$f} eq "NHR-25"){
			$ifexp=1;
			$nhr25title.="\t".$f;
			%vexp=%nhr25exp;
			$nhr25marker++;
			}
		elsif($marker{$f} eq "PHA-4"){
			$ifexp=1;
			$pha4title.="\t".$f;
			%vexp=%pha4exp;
			$pha4marker++;
			}
		open(IN,"$dir/$f") or die "cannot open file $f: $!";
		print STDERR "$f\n";
		print LOG "$f\n";
		%vblot=();
		%start=();
		%end=();
		%plusblot=();
		our %cellendflag=();
		$ct="NA";
		$thresh=$tpend{$f};
		@e2cc=();	
		@e4cc=();	
		@e8cc=();
		open(FILE,">>$cccor") or die "can not open $cccor,$!\n";
		print FILE "$f\t$tpend{$f}\t$endcell{$f}\t$ko{$f}\t";
		close(FILE);
		while(<IN>){
			chomp;
			@t=split(/,/);
			if(not exists $start{$t[1]}){
				$start{$t[1]}=$t[2];
				}
			$end{$t[1]}=$t[2];
			if($t[2]<=$thresh){
				if($t[6]>0){$plusblot{$t[1]}++;}
				if(not exists $vblot{$t[1]}){$vblot{$t[1]}=$t[1].",".$t[2].",".$t[8].",".$t[9].",".$t[10].",".$t[6];}
				else{$vblot{$t[1]}.=",".$t[1].",".$t[2].",".$t[8].",".$t[9].",".$t[10].",".$t[6];}
				}
			}
		close(IN);
		# cell cycle and expression determination
		my $outa="tmp.cc.".$$;
		my $outb="tmp.cc.bigVar.".$$;
		open(OUTA,">$outa");
		open(OUTB,">$outb");
		$ccflag=0;
		$diffc1="NA";
		$diffc2="NA";
		$p4v="NA";
		$caa="NA";
		$cap="NA";
		$caaa="NA";
		$caap="NA";
		foreach my $cell(keys %wtcc){
			#cell cycle
			if(exists $end{$cell} && $end{$cell}<=$thresh){
				$ct=($end{$cell}-$start{$cell}+1)*$cor;
				if($wtcc{$cell}<30 && $ct<=10.5){$ct="NA";} ##to eliminate editing error effect
				elsif($wtcc{$cell}>=30 && $ct<=15){$ct="NA";}
				$cellendflag{$cell}=1;
				}
				else{
					$ct="NA";
					}
			if(not exists $vcc{$cell}){$vcc{$cell}=$ct;}
			else{$vcc{$cell}=$vcc{$cell}."\t".$ct;}
			
			#checking cell cycle
			if(exists $e2flag{$cell}){push(@e2cc,$ct);}
			if(exists $e4flag{$cell}){push(@e4cc,$ct);}
			if(exists $e8flag{$cell}){push(@e8cc,$ct);}
			if($wtcc{$cell} ne "NA" && $ct ne "NA"){
				print OUTA "$cell\t$wtcc{$cell}\t$ct\n";
				$wv=$wtcc{$cell};
				$vv=$ct;
				if($vv/$wv>1.5 || $wv/$vv>1.5){print OUTB "$cell\t$wtcc{$cell}\t$ct\n";$ccflag++;}
				}
			if($cell eq 'Caa' && $ct ne "NA"){$caa=$ct;}
			if($cell eq 'Cap' && $ct ne "NA"){$cap=$ct;}
			if($cell eq 'Caaa' && $ct ne "NA"){$caaa=$ct;}
			if($cell eq 'Caap' && $ct ne "NA"){$caap=$ct;}
			if($cell eq 'P4' && $ct ne "NA"){$p4v=$ct;}
			}
		close(OUTA);
		close(OUTB);
		
		#cell expression
		if($ifexp==1){
			our @tpblots=();
			foreach $c(keys %vblot){
				if($plusblot{$c}>10){push(@tpblots,$vblot{$c});}
				}
			our %endexp=();
			lmexp(@tpblots);
			foreach my $cell(keys %wtcc){
				if(exists $endexp{$cell}){$expr=$endexp{$cell};}
				else{$expr="NA";}
				if(not exists $vexp{$cell}){$vexp{$cell}=$expr;}
				else{$vexp{$cell}.="\t".$expr;}
				}
			if($marker{$f} eq "NHR-25"){
				%nhr25exp=%vexp;
				}
			elsif($marker{$f} eq "PHA-4"){
				%pha4exp=%vexp;
				}
			}
		
		#rank tree cell cycle and calc tree cor
		$type="cc";
		$similar=rank($outa,$cccor,$f,$type,$outb,$ccflag);
		open(FILE,">>$cccor") or die "can not open $cccor,$!\n";
		$e2dif=mean(@e2cc)-$e2ave;
		$e4dif=mean(@e4cc)-$e4ave;
		$e8dif=mean(@e8cc)-$e8ave;
		print FILE "\t$e2dif\t$e4dif\t$e8dif\t";
		if($caa ne "NA" && $cap ne "NA"){$diffc1=$cap-$caa;}
		if($caaa ne "NA" && $caap ne "NA"){$diffc2=$caap-$caaa;}
		print FILE "$diffc1\t$diffc2\t$p4v\n";
		close(FILE);
		`rm $outa $outb`;
		}
	}

print CCM "$cctitle\n";
print PHA4MA "$pha4title\n";
print NHR25MA "$nhr25title\n";
foreach my $cell(keys %wtcc){
	print CCM "$cell\t$wtcc{$cell}\t$vcc{$cell}\n";
	print PHA4MA "$cell\t$wexppha4{$cell}\t$pha4exp{$cell}\n";
	print NHR25MA "$cell\t$wexpnhr25{$cell}\t$nhr25exp{$cell}\n";
	}
close(CCM);
close(PHA4MA);
close(NHR25MA);

open(REP,">$ccrep");
print REP "exp\trep1\trep2\trep.Abal\trep.ABar\trep.ABpl\trep.ABpr\trep.C\trep.E\trep.MS\trep.P3\trep.whole\tave2wt.ABal\tave2wt.ABar\tave2wt.ABpl\tave2wt.ABpr\tave2wt.C\tave2wt.E\tave2wt.MS\tave2wt.P3\tave2wt.whole\n";
$type="cc";
foreach $replicate(keys %exprep){
	@fcount=split(/;/,$exprep{$replicate});
	$reps=scalar(@fcount);
	print STDERR "$replicate\t$reps\n";
	if($reps>1){
		for(my $z=0;$z<$reps-1;$z++){
			for(my $zz=$z+1;$zz<$reps;$zz++){
				my $outa="tmp.exp.rep.".$$;
				my $outb="tmp.bigVar.rep.".$$;
				my $outd="tmp.exp.ave.".$$;
				my $oute="tmp.bigVar.ave.".$$;
				open(OUTA,">$outa");
				open(OUTB,">$outb");
				open(OUTD,">$outd");
				open(OUTE,">$oute");
				$flag1=0;
				$flag2=0;
				foreach $cell(keys %wtcc){
					@t=split(/\t/,$vcc{$cell});
					$rep1=$t[$fcount[$z]];
					$rep2=$t[$fcount[$zz]];
					$repmean=mean(($rep1,$rep2));
					if($rep1 ne "NA" && $rep2 ne "NA"){
							print OUTA "$cell\t$rep1\t$rep2\n";
							$tmp1=($rep1==0?1e-10:$rep1);
							$tmp2=($rep2==0?1e-10:$rep2);
							if($tmp1/$tmp2>1.5 || $tmp2/$tmp1>1.5){print OUTB "$cell\t$rep1\t$rep2\n";$flag1++;}
							}
					if($wtcc{$cell} ne "NA" && $repmean!=0){
						print OUTD "$cell\t$wtcc{$cell}\t$repmean\n";
						$wv=($wtcc{$cell}==0?1e-10:$wtcc{$cell});
						$vv=($repmean==0?1e-10:$repmean);
						if($vv/$wv>1.5 || $wv/$vv>1.5){print OUTE "$cell\t$wtcc{$cell}\t$repmean\n";$flag2++;}
						}
					}
				close(OUTA);
				close(OUTB);
				close(OUTD);
				close(OUTE);
				$cors=repcor($outa,$outd,$outb,$oute,$flag1,$flag2,$replicate,$exprepname{$fcount[$z]},$exprepname{$fcount[$z+1]},$type);
				print REP "$replicate\t$exprepname{$fcount[$z]}\t$exprepname{$fcount[$zz]}\t$cors\n";
				`rm $outa $outb $outd $oute`;
				}
			}
		}
	}
close(REP);

$type="exp";
#pha4 expression
if($pha4marker>0){
	$nme=normexp($mat1);
	our %exps=();
	open(IN,"$nme") or die "can not open $nme,$!\n";
	print PHA4MB "$pha4title\n";
	while(<IN>){
		chomp;
		@t=split(/\t/);
		$id=$t[0];
		shift(@t);
		$exps{$id}=join("\t",@t);
		print PHA4MB "$id\t$wexppha4{$id}\t$exps{$id}\n";
		}
	close(PHA4MB);
	close(IN);
	`rm $nme`;
	($t1,$t2,@fs)=split(/\t/,$pha4title);
	$fc1=scalar(@fs);
	%wtexp=%wexppha4;
	%exprep=();
	%exprepname=();
	for($i=0;$i<$fc1;$i++){
		open(EXP,">>$expcor") or die "can not open $expcor,$!\n";
		print EXP "\n$fs[$i]\t$tpend{$fs[$i]}\t$endcell{$fs[$i]}\t$ko{$fs[$i]}\t";
		close(EXP);
		
		$expname=substr($fs[$i],0,-7);
		$expname=~s/^CD//;
		if(not exists $exprep{$expname}){$exprep{$expname}=$i;}
		else{$exprep{$expname}.=";".$i;}
		$exprepname{$i}=$fs[$i];

		my $outa="tmp.exp.".$$;
		my $outb="tmp.exp.bigVar.".$$;
		open(OUTA,">$outa");
		open(OUTB,">$outb");
		$expflag=0;
		@lost=();
		@raised=();
		#checking cell expression
		foreach my $cell(keys %exps){
			@t=split(/\t/,$exps{$cell});
			$expr=$t[$i];
			if($wtexp{$cell} ne "NA" && $expr ne "NA"){
				print OUTA "$cell\t$wtexp{$cell}\t$expr\n";
				$wv=$wtexp{$cell};
				$vv=$expr;
				if($wv-$vv>0.15 && $wv>1 && $vv<1){
					print OUTB "$cell\t$wtexp{$cell}\t$expr\n";
					$expflag++;
					push(@lost,$cell);
					}
				elsif($vv-$wv>0.15 && $vv>1 && $wv<1){
					print OUTB "$cell\t$wtexp{$cell}\t$expr\n";
					$expflag++;
					push(@raised,$cell);
					}
				}
			}
		close(OUTA);
		close(OUTB);
		# $similar=rank($outa,$cccor,$f,$type,$outb,$ccflag);
		rank($outa,$expcor,$fs[$i],$type,$outb,$expflag);
		# plot bigVar cell
		if($expflag<3){$expflag=0;}
		if($expflag>0){
			print BV "$fs[$i]\t";
			$a=join(", ",@raised);
			print BV "$a\t";
			$a=join(", ",@lost);
			print BV "$a\n";	
			}
		`rm $outa $outb`;
		# `mv $outa $fs[$i].outa`;
		# `mv $outb $fs[$i].outb`;
		}

	open(REP,">$exprrep");
	print REP "exp\trep1\trep2\trep.Abal\trep.ABar\trep.ABpl\trep.ABpr\trep.C\trep.E\trep.MS\trep.P3\trep.whole\tave2wt.ABal\tave2wt.ABar\tave2wt.ABpl\tave2wt.ABpr\tave2wt.C\tave2wt.E\tave2wt.MS\tave2wt.P3\tave2wt.whole\n";
	foreach $replicate(keys %exprep){
		@fcount=split(/;/,$exprep{$replicate});
		$reps=scalar(@fcount);
		print STDERR "$replicate\t$reps\n";
		if($reps>1){
			for(my $z=0;$z<$reps-1;$z++){
				for(my $zz=$z+1;$zz<$reps;$zz++){
					my $outa="tmp.exp.rep.".$$;
					my $outb="tmp.bigVar.rep.".$$;
					my $outd="tmp.exp.ave.".$$;
					my $oute="tmp.bigVar.ave.".$$;
					open(OUTA,">$outa");
					open(OUTB,">$outb");
					open(OUTD,">$outd");
					open(OUTE,">$oute");
					$flag1=0;
					$flag2=0;
					foreach $cell(keys %wtexp){
						@t=split(/\t/,$exps{$cell});
						$rep1=$t[$fcount[$z]];
						$rep2=$t[$fcount[$zz]];
						$repmean=mean(($rep1,$rep2));
						if($rep1 ne "NA" && $rep2 ne "NA"){
								print OUTA "$cell\t$rep1\t$rep2\n";
								$tmp1=$rep1;
								$tmp2=$rep2;
								if(($tmp2-$tmp1>0.15 && $tmp1<1 && $tmp2>1) || ($tmp1-$tmp2>0.15 && $tmp1>1 && $tmp2<1)){
									print OUTB "$cell\t$rep1\t$rep2\n";
									$flag1++;
									}
								}
						if($wtexp{$cell} ne "NA" && $repmean!=0){
							print OUTD "$cell\t$wtexp{$cell}\t$repmean\n";
							$wv=$wtexp{$cell};
							$vv=$repmean;
							if($wv-$vv>0.15 && $wtexp{$cell}>1){
								print OUTE "$cell\t$wtexp{$cell}\t$repmean\n";
								$flag2++;
								}
							elsif($vv-$wv>0.15 && $repmean>1){
								print OUTE "$cell\t$wtexp{$cell}\t$repmean\n";
								$flag2++;
								}
							}
						}
					close(OUTA);
					close(OUTB);
					close(OUTD);
					close(OUTE);
					$cors=repcor($outa,$outd,$outb,$oute,$flag1,$flag2,$replicate,$exprepname{$fcount[$z]},$exprepname{$fcount[$z+1]},$type);
					print REP "$replicate\t$exprepname{$fcount[$z]}\t$exprepname{$fcount[$zz]}\t$cors\n";
					`rm $outa $outb $outd $oute`;
					}
				}
			}
		}
	}

#nhr-25 expression
if($nhr25marker>1){
	$nme=normexp($mat3);
	our %exps=();
	open(IN,"$nme") or die "can not open $nme,$!\n";
	print NHR25MB "$nhr25title\n";
	while(<IN>){
		chomp;
		@t=split(/\t/);
		$id=$t[0];
		shift(@t);
		$exps{$id}=join("\t",@t);
		print NHR25MB "$id\t$wexpnhr25{$id}\t$exps{$id}\n";
		}
	close(NHR25MB);
	close(IN);
	`rm $nme`;
	($t1,$t2,@fs)=split(/\t/,$nhr25title);
	$fc2=scalar(@fs);
	%wtexp=%wexpnhr25;
	%exprep=();
	%exprepname=();
	for($i=0;$i<$fc2;$i++){
		open(EXP,">>$expcor") or die "can not open $expcor,$!\n";
		print EXP "\n$fs[$i]\t$tpend{$fs[$i]}\t$endcell{$fs[$i]}\t$ko{$fs[$i]}\t";
		close(EXP);
		
		$expname=substr($fs[$i],0,-7);
		$expname=~s/^CD//;
		if(not exists $exprep{$expname}){$exprep{$expname}=$i;}
		else{$exprep{$expname}.=";".$i;}
		$exprepname{$i}=$fs[$i];
		
		my $outa="tmp.exp.".$$;
		my $outb="tmp.exp.bigVar.".$$;
		open(OUTA,">$outa");
		open(OUTB,">$outb");
		$expflag=0;
		@lost=();
		@raised=();
		#checking cell expression
		foreach my $cell(keys %exps){
			@t=split(/\t/,$exps{$cell});
			$expr=$t[$i];
			if($wtexp{$cell} ne "NA" && $expr ne "NA"){
				print OUTA "$cell\t$wtexp{$cell}\t$expr\n";
				$wv=$wtexp{$cell};
				$vv=$expr;
				if($wv-$vv>0.15 && $wv>1 && $vv<1){
					print OUTB "$cell\t$wtexp{$cell}\t$expr\n";
					$expflag++;
					push(@lost,$cell);
					}
				elsif($vv-$wv>0.15 && $vv>1 && $wv<1){
					print OUTB "$cell\t$wtexp{$cell}\t$expr\n";
					$expflag++;
					push(@raised,$cell);
					}
				}
			}
		close(OUTA);
		close(OUTB);
		rank($outa,$expcor,$fs[$i],$type,$outb,$expflag);
		# plot bigVar cell
		if($expflag<3){$expflag=0;}
		if($expflag>0){
			print BV "$fs[$i]\t";
			$a=join(", ",@raised);
			print BV "$a\t";
			$a=join(", ",@lost);
			print BV "$a\n";	
			}
		`rm $outa $outb`;
		}
	close(BV);

	foreach $replicate(keys %exprep){
		@fcount=split(/;/,$exprep{$replicate});
		$reps=scalar(@fcount);
		print STDERR "$replicate\t$reps\n";
		if($reps>1){
			for(my $z=0;$z<$reps-1;$z++){
				for(my $zz=$z+1;$zz<$reps;$zz++){
					my $outa="tmp.exp.rep.".$$;
					my $outb="tmp.bigVar.rep.".$$;
					my $outd="tmp.exp.ave.".$$;
					my $oute="tmp.bigVar.ave.".$$;
					open(OUTA,">$outa");
					open(OUTB,">$outb");
					open(OUTD,">$outd");
					open(OUTE,">$oute");
					$flag1=0;
					$flag2=0;
					foreach $cell(keys %wtexp){
						@t=split(/\t/,$exps{$cell});
						$rep1=$t[$fcount[$z]];
						$rep2=$t[$fcount[$zz]];
						$repmean=mean(($rep1,$rep2));
						if($rep1 ne "NA" && $rep2 ne "NA"){
								print OUTA "$cell\t$rep1\t$rep2\n";
								$tmp1=$rep1;
								$tmp2=$rep2;
								if(($tmp2-$tmp1>0.2 && $tmp1<1 && $tmp2>1) || ($tmp1-$tmp2>0.2 && $tmp1>1 && $tmp2<1)){
									print OUTB "$cell\t$rep1\t$rep2\n";
									$flag1++;
									}
								}
						if($wtexp{$cell} ne "NA" && $repmean!=0){
							print OUTD "$cell\t$wtexp{$cell}\t$repmean\n";
							$wv=$wtexp{$cell};
							$vv=$repmean;
							if($wv-$vv>0.2 && $wtexp{$cell}>1){
								print OUTE "$cell\t$wtexp{$cell}\t$repmean\n";
								$flag2++;
								}
							elsif($vv-$wv>0.2 && $repmean>1){
								print OUTE "$cell\t$wtexp{$cell}\t$repmean\n";
								$flag2++;
								}
							}
						}
					close(OUTA);
					close(OUTB);
					close(OUTD);
					close(OUTE);
					$cors=repcor($outa,$outd,$outb,$oute,$flag1,$flag2,$replicate,$exprepname{$fcount[$z]},$exprepname{$fcount[$z+1]},$type);
					print REP "$replicate\t$exprepname{$fcount[$z]}\t$exprepname{$fcount[$zz]}\t$cors\n";
					`rm $outa $outb $outd $oute`;
					}
				}
			}
		}
	}
close(REP);

close(LOG);

sub rank(){
	my ($infile,$outfile,$cdfile,$datatype,$varfile,$varflag)=@_;
	open(FILE,">>$outfile");
	my %small=();
	my %big=();
	my $rankbig="";
	my $ranksmall="";
	my $stablebig="";
	my $stablesmall="";
	my $corscore="";
	$small{'ABal'}=meanCanberra($infile,@abaltree);
	$small{'ABar'}=meanCanberra($infile,@abartree);
	$small{'ABpl'}=meanCanberra($infile,@abpltree);
	$small{'ABpr'}=meanCanberra($infile,@abprtree);
	$big{'C'}=meanCanberra($infile,@ctree);
	$small{'Ca'}=meanCanberra($infile,@catree);
	$small{'Cp'}=meanCanberra($infile,@cptree);
	$big{'E'}=meanCanberra($infile,@etree);
	$small{'Ea'}=meanCanberra($infile,@eatree);
	$small{'Ep'}=meanCanberra($infile,@eptree);
	$big{'MS'}=meanCanberra($infile,@mstree);
	$small{'MSa'}=meanCanberra($infile,@msatree);
	$small{'MSp'}=meanCanberra($infile,@msptree);
	$small{'D'}=meanCanberra($infile,@dtree);
	$big{'P3'}=meanCanberra($infile,@p3tree);
	$rsqall=meanCanberra($infile,@all);
	my @bigt=($small{'ABal'},$small{'ABar'},$small{'ABpl'},$small{'ABpr'},$big{'C'},$big{'E'},$big{'MS'});
	my @smallt=($small{'ABal'},$small{'ABar'},$small{'ABpl'},$small{'ABpr'},$small{'Ca'},$small{'Cp'},$small{'Ea'},$small{'Ep'},$small{'MSa'},$small{'MSp'});
	my %namebig=('0'=>'ABal','1'=>'ABar','2'=>'ABpl','3'=>'ABpr',
				'4'=>'C','5'=>'E','6'=>'MS','7'=>'P3');
	my %namesmall=('0'=>'ABal','1'=>'ABar','2'=>'ABpl','3'=>'ABpr',
				'4'=>'Ca','5'=>'Cp','6'=>'Ea','7'=>'Ep','8'=>'MSa','9'=>'MSp','10'=>'D');
	$corscore=join("\t",@smallt);
	$corscore.="\t".$small{'D'};
	$exp=$cdfile;
	$exp=~s/^CD//;
	$exp=~s/\.csv$//;
	my $size=scalar(@bigt);
	for(my $i=0;$i<$size;$i++){
		@temparray=@bigt;
		$tempv=splice(@temparray,$i,1);
		$mint=min @temparray;
		$maxt=max @temparray;
		if($mint-$tempv>0.3 ||($mint>0.9 && $tempv<0.8)){
			$rankbig=$namebig{$i};
			$ifsep=plottree($infile,$rsqall,$tempv,$exp,$namebig{$i},$datatype);
			if($ifsep==1){$rankbig.=" (distinct)";}
			}
		if($tempv-$maxt>0.3 ||($tempv>0.9 && $maxt<0.8)){
			$stablebig=$namebig{$i};
			$ifsep=plottree($infile,$rsqall,$tempv,$exp,$namebig{$i},$datatype);
			if($ifsep==1){$stablebig.=" (distinct)";}
			}
		}
	my $size=scalar(@smallt);
	for(my $i=0;$i<$size;$i++){
		@temparray=@smallt;
		$tempv=splice(@temparray,$i,1);
		$mint=min @temparray;
		$maxt=max @temparray;
		if($mint-$tempv>0.3 ||($mint>0.9 && $tempv<0.8)){
			$ranksmall=$namesmall{$i};
			$ifsep=plottree($infile,$rsqall,$tempv,$exp,$namesmall{$i},$datatype);
			if($ifsep==1){$ranksmall.=" (distinct)";}
			}
		if($tempv-$maxt>0.3 ||($tempv>0.9 && $maxt<0.8)){
			$stablesmall=$namesmall{$i};
			$ifsep=plottree($infile,$rsqall,$tempv,$exp,$namesmall{$i},$datatype);
			if($ifsep==1){$stablesmall.=" (distinct)";}
			}		
		}
	my @et=($big{'E'},$small{'Ea'},$small{'Ep'});
	my @mst=($big{'MS'},$small{'MSa'},$small{'MSp'});
	my @abt=($small{'ABal'},$small{'ABar'},$small{'ABpl'},$small{'ABpr'});
	my @ct=($big{'C'},$small{'Ca'},$small{'Cp'});
	$comems="";
	$comcms="";
	$comcab="";
	$comab="";
	if(min(@et)-max(@mst)>0.2){$comems="MS changed";}
	if(min(@mst)-max(@et)>0.2){$comems="E changed";}		
	if(min(@ct)-max(@mst)>0.2){$comcms="MS changed";}
	if(min(@mst)-max(@ct)>0.2){$comcms="C changed";}
	if(min(@ct)-max(@abt)>0.2){$comcab="AB changed";}
	if(min(@abt)-max(@ct)>0.2){$comcab="C changed";}
	my %nameab=('0'=>'ABal','1'=>'ABar','2'=>'ABpl','3'=>'ABpr');
	$size=scalar(@abt);
	for(my $i=0;$i<$size;$i++){
		@temparray=@smallt;
		$tempv=splice(@temparray,$i,1);
		$mint=min @temparray;
		if($mint-$tempv>0.2){
			$comab=$nameab{$i}." changed";
			}
		}
	print FILE "$small{'ABal'}\t$small{'ABar'}\t$small{'ABpl'}\t$small{'ABpr'}\t";
	print FILE "$big{'C'}\t$small{'Ca'}\t$small{'Cp'}\t";
	print FILE "$big{'E'}\t$small{'Ea'}\t$small{'Ep'}\t";
	print FILE "$big{'MS'}\t$small{'MSa'}\t$small{'MSp'}\t";
	print FILE "$small{'D'}\t$big{'P3'}\t$rsqall\t";
	print FILE "$ranksmall\t$rankbig\t";
	print FILE "$stablesmall\t$stablebig\t";
	print FILE "$comems\t$comcms\t$comcab\t$comab";
	close(FILE);
	# plot bigVar cell
	if($varflag>0){
		plotbigvar($infile,$varfile,$rsqall,$exp,$datatype);
		}
	return($corscore);
	}

sub repcor(){
	my ($infile1,$infile2,$varfile1,$varfile2,$varflag1,$varflag2,$name,$name1,$name2,$datatype)=@_;
	@rv=();
	$rv[0]=meanCanberra($infile1,@abaltree);
	$rv[1]=meanCanberra($infile1,@abartree);
	$rv[2]=meanCanberra($infile1,@abpltree);
	$rv[3]=meanCanberra($infile1,@abprtree);
	$rv[4]=meanCanberra($infile1,@ctree);
	$rv[5]=meanCanberra($infile1,@etree);
	$rv[6]=meanCanberra($infile1,@mstree);
	$rv[7]=meanCanberra($infile1,@p3tree);
	$rv[8]=meanCanberra($infile1,@all);
	$rv[9]=meanCanberra($infile2,@abaltree);
	$rv[10]=meanCanberra($infile2,@abartree);
	$rv[11]=meanCanberra($infile2,@abpltree);
	$rv[12]=meanCanberra($infile2,@abprtree);
	$rv[13]=meanCanberra($infile2,@ctree);
	$rv[14]=meanCanberra($infile2,@etree);
	$rv[15]=meanCanberra($infile2,@mstree);
	$rv[16]=meanCanberra($infile2,@p3tree);
	$rv[17]=meanCanberra($infile2,@all);
	$rvs=join("\t",@rv);
	if($varflag1>0){
		plotrep($infile1,$varfile1,$rv[8],$name,$name1,$name2,$datatype);
		}
	if($varflag2>0){
		plotbigvar($infile2,$varfile2,$rv[17],$name,$datatype);
		}
	return($rvs);
	}
	
sub plottree(){
	my ($infile,$rall,$rtree,$expn,$treen,$pren)=@_;
	my @array=split(/;/,$treearray{$treen});
	my $intree="tmp.tree.".$$;
	my $rf="tmp.".$$.".R";
	my $rval=0;
	my $inflag=0;
	my %tmphash=();
	foreach $tmpv(@array){
		$tmphash{$tmpv}=1;
		}
	open(IN,$infile);
	open(A,">$intree");
	my @treed=();
	my @otherd=();
	while(<IN>){
		chomp;
		my @a=split(/\t/);
		if(exists $tmphash{$a[0]}){
			print A "$_\n";
			$inflag++;
			push(@treed,abs($a[1]-$a[2]));
			}
		else{
			push(@otherd,abs($a[1]-$a[2]));
			}
		}
	close(IN);
	close(A);
	if($inflag>5){
		# if(($tmpave2-$tmpave1>6.5 && $maxtreedif>$maxdif && $mintreedif>$mindif) or ($tmpave1-$tmpave2>6.5 && $maxtreedif<$maxdif && $mintreedif<$mindif)){
		my $tmpave1=mean(@treed);
		my $tmpave2=mean(@otherd);
		my $tmpmd1=median(@treed);
		my $tmpmd2=median(@otherd);
		print LOG "$treen\t$inflag\t$tmpave2\t$tmpave1\t$tmpmd2\t$tmpmd1\n";
		if(($tmpave2-$tmpave1>4.5 && $tmpmd2-$tmpmd1>3) or ($tmpave1-$tmpave2>4.5 && $tmpmd1-$tmpmd2>3)){
			$rval=1;
			}
		my $outp=$pren.".tree.".$treen.".".$expn.".pdf";
		my $outpa=$pren.".tree.".$treen.".".$expn.".png";
		my $path=`pwd`;
		chomp $path;
		open(A,">$rf");
		print A qq(setwd("$path")
read.table("$infile",sep="\t")->x
read.table("$intree",sep="\t")->xx
reg1<-lm(x[,2]~x[,3])
reg2<-lm(xx[,2]~xx[,3])
flag<-$rval
axismin<-min(min(x[,2],na.rm=TRUE),min(x[,3],na.rm=TRUE))
axismax<-max(max(x[,2],na.rm=TRUE),max(x[,3],na.rm=TRUE))
rsq<-paste("Similarity=$rall\nSimilarity.$treen=$rtree")
if(flag==1){
	name<-paste("dist.","$outp",sep="")
	namea<-paste("dist.","$outpa",sep="")
	}
if(flag==0){
	name<-paste("$outp")
	namea<-paste("$outpa")
	}
pdf(file=name)
plot(x[,3],x[,2],col="black",pch=20,ylab="Wild Type",xlab="$expn",xlim=c(axismin,axismax),ylim=c(axismin,axismax))
points(xx[,3],xx[,2],col="red",pch=20)
text(xx[,3],xx[,2],labels=as.character(xx[,1]),pos=4,col="red")
legend("topleft",rsq,bty = "n")
abline(reg1,col="black")
abline(reg2,col="red")
dev.off()
png(file=namea,width=1000,height=1000)
plot(x[,3],x[,2],col="black",pch=20,ylab="Wild Type",xlab="$expn",xlim=c(axismin,axismax),ylim=c(axismin,axismax))
points(xx[,3],xx[,2],col="red",pch=20)
text(xx[,3],xx[,2],labels=as.character(xx[,1]),pos=4,col="red")
legend("topleft",rsq,bty = "n")
abline(reg1,col="black")
abline(reg2,col="red")
dev.off()
);
		close(A);
		`Rscript $rf`;
		`rm $intree $rf`;
		}
	return($rval);
	}
	
sub plotbigvar(){
	my ($infile,$varfile,$rall,$expn,$pren)=@_;
	my $rf="tmp.".$$.".R";
	my $outp=$pren.".bigVar.".$expn.".png";
	my $path=`pwd`;
	chomp $path;
	open(OUTA,">$rf");
	print OUTA qq(setwd("$path")
read.table("$infile",sep="\t")->x
read.table("$varfile",sep="\t")->xx
axismin<-min(min(x[,2],na.rm=TRUE),min(x[,3],na.rm=TRUE))
axismax<-max(max(x[,2],na.rm=TRUE),max(x[,3],na.rm=TRUE))
rsq<-paste("Similarity=$rall")
name<-paste("$outp")
png(file=name,width=1000,height=1000)
reg1<-lm(x[,2]~x[,3])
plot(x[,3],x[,2],col="black",pch=20,ylab="Wild Type",xlab="$expn",xlim=c(axismin,axismax),ylim=c(axismin,axismax))
points(xx[,3],xx[,2],col="red",pch=20)
text(xx[,3],xx[,2],labels=as.character(xx[,1]),pos=4,col="red")
legend("topleft",rsq,bty = "n")
abline(reg1,col="black")
dev.off()
);
	close(OUTA);
	`Rscript $rf`;
	`rm $rf`;
	}

sub mean{
	my(@data)=@_;
	if (not @data) {
		die("Empty array\n");
	}
	my $total=0;
	my $count=0;
	foreach (@data) {
		if($_ ne "NA"){
			$total+=$_;
			$count++;
			}
	}
	if($total==0 && $count==0){$count=1;}
	my $average=$total/$count;
	return $average;
}

sub plotrep(){
	my ($cctfile,$varfile,$rall,$expn,$exp1,$exp2,$pren)=@_;
	my $rf="tmp.".$$.".R";
	my $outp=$pren.".rep.bigVar.".$expn.".png";
	my $path=`pwd`;
	chomp $path;
	open(A,">$rf");
	print A qq(setwd("$path")
read.table("$cctfile",sep="\t")->x
read.table("$varfile",sep="\t")->xx
axismin<-min(min(x[,2],na.rm=TRUE),min(x[,3],na.rm=TRUE))
axismax<-max(max(x[,2],na.rm=TRUE),max(x[,3],na.rm=TRUE))
rsq<-paste("Similarity=$rall")
name<-paste("$outp")
png(file=name,width=1000,height=1000)
reg1<-lm(x[,2]~x[,3])
plot(x[,3],x[,2],col="black",pch=20,main="$expn",xlab="$exp2",ylab="$exp1",xlim=c(axismin,axismax),ylim=c(axismin,axismax))
points(xx[,3],xx[,2],col="red",pch=20)
text(xx[,3],xx[,2],labels=as.character(xx[,1]),pos=4,col="red")
legend("topleft",rsq,bty = "n")
abline(reg1,col="black")
dev.off()
);
	close(A);
	`Rscript $rf`;
	`rm $rf`;
	}

sub lmexp(){
	my $rf="tmp.".$$.".R";
	my $infile="tmp.cell.exp.".$$;
	my $endfile="tmp.cell.end.".$$;
	my $rout="tmp.out.".$$;
	open(A,">$infile");
	open(B,">$endfile");
	my $rowflag=1;
	my $rowend=1;
	foreach $exps(@_){
		my @ax=split(/,/,$exps);
		my $size1=scalar(@ax);
		my $size=$size1/6;
		my $sum=0;
		my $cellname=$ax[0];
		if($size>10 && $plusblot{$ax[0]}>8){
			for(my $k=3;$k<$size-2;$k++){
				my $ka=($k-1)*6+5;
				my $kb=$k*6+5;
				$sum+=abs($ax[$kb]-$ax[$ka]);
				}
			my $avediff=$sum/($size-5);
			if($avediff<200){$avediff=$avediff*2;}
			elsif($avediff<500 && $avediff>=200){$avediff=$avediff*1.5;}
			my $last=($ax[($size-3)*6+5]+$ax[($size-4)*6+5]+$ax[($size-5)*6+5])/3;
			my $nacount=0;
			for($k=4;$k<$size;$k++){
				$ka=$k*6-1;
				$kb=$k*6+5;
				while($ax[$ka] eq "NA"){
					$ka=$ka-6;
					if($ax[$ka] ne "NA"){last;}
					}
				if(($ax[$kb]-$ax[$ka])/2>$avediff && $ax[$kb]/2>$last && $ax[$kb]>2000){print LOG "$cellname\t$avediff\t$last\t$ax[$ka]\t$ax[$kb]\n";$ax[$kb]="NA";$nacount++;}
				elsif($ax[$kb]/3>$ax[$ka] && $ax[$ka]>500 && $ax[$kb]/2>$last && $ax[$kb]>2000){print LOG "$cellname\t$avediff\t$last\t$ax[$ka]\t$ax[$kb]\n";$ax[$kb]="NA";$nacount++;}
				}
			my @xx=();
			if($size-$nacount>6){
				for($k=0;$k<$size;$k++){
					push(@xx,$ax[$k*6+1],$ax[$k*6+2],$ax[$k*6+3],$ax[$k*6+4],$ax[$k*6+5]);
					}
				}
			print B "$cellname\t$rowflag\t";
			for($k=2;$k<$size-2;$k++){
				$rowend=$rowflag;
				print A "$xx[$k*5]\t$xx[$k*5+1]\t$xx[$k*5+2]\t$xx[$k*5+3]\t$xx[$k*5+4]\n";
				$rowflag++;
				}
			if(not exists $cellendflag{$cell}){
				$rowend=$rowflag;
				print A "$xx[$k*5]\t$xx[$k*5+1]\t$xx[$k*5+2]\t$xx[$k*5+3]\t$xx[$k*5+4]\n";
				$rowflag++;
				$k++;
				$rowend=$rowflag;
				print A "$xx[$k*5]\t$xx[$k*5+1]\t$xx[$k*5+2]\t$xx[$k*5+3]\t$xx[$k*5+4]\n";
				$rowflag++;
				print B "$rowend\t$xx[$k*5]\t$xx[$k*5+1]\t$xx[$k*5+2]\t$xx[$k*5+3]\n";
				}
			else{
				print B "$rowend\t$xx[($size-1)*5]\t$xx[($size-1)*5+1]\t$xx[($size-1)*5+2]\t$xx[($size-1)*5+3]\n";
				}
			}
		}
		close(A);
		close(B);
		my $path=`pwd`;
		chomp $path;
		open(A,">$rf");
		print A qq(setwd("$path")
library(glmnet)
read.table("$infile",sep="\t")->x
read.table("$endfile",sep="\t")->y
nrow(y)->nr
score<-rep("NA",nr)
for (i in 1:nr){
	na.omit(x[y[i,2]:y[i,3],1:5])->xx
	lasso_fit<-glmnet(as.matrix(xx[,1:4]),as.matrix(xx[,5]),alpha=1)
	score[i]<-predict(lasso_fit,as.matrix(y[i,4:7]),s=1)
	}
write.table(data.frame(y[,1],score),file="$rout",sep="\t",quote=F)
);
		close(A);
		`Rscript $rf`;
		open(A,$rout);
		<A>;
		while(<A>){
			chomp;
			@t=split(/\t/);
			if($t[2]<=10){$t[2]="NA";}
			$endexp{$t[1]}=$t[2];
			}
		close(A);
		`rm $rf $infile $endfile $rout`;
	}
		
sub normexp(){
	my ($f)=@_;
	my $rf="tmp.".$$.".R";
	my $rout="tmp.R.out.".$$;
	my $path=`pwd`;
	chomp $path;
	open(A,">$rf");
	print A qq(setwd("$path")
library(limma)
library(vsn)
as.matrix(read.table("$f",sep="\t",row.names=1,header=T))->x
normalizeVSN(x[,-1])->xx
ncol(xx)->max
xx->y
for(i in 1:max){
	xx[,i]/quantile(na.omit(xx[,i]))[4]->y[,i]
	}
write.table(y,file="$rout",sep="\t",quote=F,row.names=T,col.names=F)
);
	close(A);
	`Rscript $rf`;
	`rm $rf`;
	return($rout);
	}

sub meanCanberra(){
	my ($infile,@array)=@_;
	my %tmphash=();
	foreach $tmpv(@array){
		$tmphash{$tmpv}=1;
		}
	my $tmpcount=0;
	my $tmpsum=0;
	my $tmpmean=0;
	open(IN,$infile);
	while(<IN>){
		chomp;
		my @a=split(/\t/);
		if(exists $tmphash{$a[0]}){
			$tmpsum+=abs(($a[1]-$a[2])/($a[1]+$a[2]));
			$tmpcount++;
			}
		}
	close(IN);
	if($tmpcount==0){$tmpcount=1;}
	$tmpmean=$tmpsum/$tmpcount;
	my $tmpsimilarity=(1-$tmpmean*2)**2;
	return($tmpsimilarity);
	}
