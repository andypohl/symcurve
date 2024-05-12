#!/usr/bin/perl 

#Program SymCurv_prediction.pl
#Version_Y: 23/03/2009
# ULTIMATE PROGRAM (SymCurv_v2.0)
# 1. Calculates Curvature from a fasta file, then stores it in an array and goes on with Symcurv calculations
# 2. Calculates Symmetry of Curvature around local minima. Score is calculated as 1/SUM(ABS(CURV(I-M)-CURV(I+M)))
# 3. Makes overlapping nucleosomes calls based on the SymCurv values
# 4. Deduces non-overlapping nucleosomes with the use of a greedy algorithm
# 5. Produces distinct outputs including
#	a. Curvature
#	b. SymCurv
#	c. Overlapping Calls in gff format
#	d. Non-overlapping calls in gff

$start = time (); #starting timer;

#Creating a directory unless it already exists
mkdir ("SYMCURV_OUT", 0777);
#
# Opening of output files
# Files used
open OUT0, ">SYMCURV_OUT/out_curv.dat";
open OUT1, ">SYMCURV_OUT/out_symcurv.dat";
open OUT11, ">SYMCURV_OUT/out_symcurv_aggr.dat";
open OUT12, ">SYMCURV_OUT/out_symcurv_avr.dat";
open OUT2, ">SYMCURV_OUT/out_init_calls_nuc.gff";
open OUT3, ">SYMCURV_OUT/out_init_calls_dnase.gff";
open OUT4, ">SYMCURV_OUT/out_fin_calls_nuc.gff";
open OUT5, ">SYMCURV_OUT/out_fin_calls_dnase.gff";

#
#Initial Input Settings
$title = ''; # title of the sequence under analysis. Is carried in the FASTA file title
$offset = 0; # initial coordinate of the sequence. Is carried in the FASTA file title
$strand = ''; # strand of the sequence under analyis. If "-" the revcomp is set as sequence
$count = 0; # number of sequences
$sequence = ''; # sequence as string
@DNA = ''; # sequence as array
@NumDNA = ''; # numerical array of DNA sequence (A=0 , T=1 , G=2 , C=3)
#

#
#Reading sequence
# Reading of file: Obtaining sequence
	while (<ARGV>) {
	#Step 1. Start of sequence with header line ">..."
		if (/>(.+):(.+):([\+-])/){
		$sequence = ''; #initializing sequence
		$title = $1;
		$offset = $2;
		$strand = $3;
		}
		if ((/>(.+)/) and not (/:/)) { # if file has not standard title line
		$sequence='';
		$title = $1;
		$offqset = 1;
		$strand = "+";
		}  
	#Step 2. Retrieval of sequence
               	if (($_ !~ />/) and ($_ !~ /#/)){
		chomp;               	
		$sequence .= $_;
     		}
	#Step 3. End of sequence denoted by #
		if (/#/) {
		#if (length($sequence) >= 150) {		
		$count++; # counter of sequences incremented
		#Calling Reading Sequence Subroutine		
		@DNA = (); # sequence as array
		@NumDNA= (); #sequence as numerical array
#
READFASTA($sequence);
#print "number of bases ",scalar(@DNA),"\n";
#print "****************\n";

# Calling curvature subroutine
CURVATURE(@NumDNA);

# printing curvature array
for ($i=1; $i<=scalar(@Curv_nuc); $i++) {
if (defined $Curv_nuc[$i]) {
	print OUT0 $title,"\t",$i+$offset,"\t",$DNA[$i],"\t",$Curv_nuc[$i],"\t",$Curv_simple[$i],"\n";
	}
}
#print " -> Values of predicted curvature for ",scalar(@Curv_nuc)," nucleotides stored in: out_curv.dat\n";

# Calling SymCurv subroutine for both sets of parameters

# 1. "stationary" state
SYMCURV(@Curv_nuc);
@symcurv_nuc=@symcurv;
$average_nuc=$average;

# 2. "activated" state
SYMCURV(@Curv_simple);
@symcurv_dnase=@symcurv;
$average_dnase=$average;

# printing SymCurv array 
for ($i=1; $i<=scalar(@symcurv_nuc); $i++) {
if (defined $symcurv_nuc[$i]) {
	print OUT1 $title,"\t",$i+$offset,"\t",$symcurv_nuc[$i],"\t",$symcurv_dnase[$i],"\n";
	}
}
#print " -> Values of Symmetry of curvature for ",scalar(@symcurv)," nucleotides stored in: out_symcurv.dat\n";
print  OUT11 $title,":",$offset,":",$strand,":",scalar(@Curv_nuc),"\t",$average_nuc,"\t",$average_dnase,"\n"; # printing aggragates of SymCurv for both parameters

#NEW# Version v22 FEATURE
#Calculating a running average for SymCurv to produce non-zero values 

#window size for which the averaging is done
$avwindow = 147; #size of a nucleosome 
for ($i=1; $i<=scalar(@symcurv_nuc); $i++) {
#if (defined $symcurv_nuc[$i]) {
		for ($j=$i-73; $j<=$i+73;$j++) {
		$symcurv_avnuc[$i]+=$symcurv_nuc[$j];
		$symcurv_avdnase[$i]+=$symcurv_dnase[$j];
		}
#	}
	print OUT12 $title,"\t",$i+$offset,"\t",$symcurv_avnuc[$i],"\t",$symcurv_avdnase[$i],"\n";
}


#printing overlapping nucleosome calls. At the same time overlapping arrays are created to be used in the next step

# 1. "stationary" state
my %overlapping_calls_nuc = ();
for ($dyad=1; $dyad<=scalar(@symcurv_nuc); $dyad++) {
	if ((defined $symcurv_nuc[$dyad]) and ($symcurv_nuc[$dyad]>0) and ($symcurv_nuc[$dyad]< 100) and ($dyad-73>0) and ($dyad+73 < length($sequence))) {
	print OUT2 $title,"\tevidence\t","stat_nucleosome\t",$dyad-73+$offset,"\t",$dyad+73+$offset,"\t",$symcurv_nuc[$dyad],"\t+\t.\t",join('', @DNA[$dyad-73..$dyad+73]),"\n";
	$overlapping_calls_nuc{$dyad}=$symcurv_nuc[$dyad];  #creating hash non overlapping calls
	}
	if ((defined $symcurv_nuc[$dyad]) and ($symcurv_nuc[$dyad]>=100) and ($dyad-73>0) and ($dyad+73 < length($sequence))) {
	print OUT2 $title,"\tevidence\t","stat_nucleosome\t",$dyad-73+$offset,"\t",$dyad+73+$offset,"\t100.0\t+\t.\t",join('', @DNA[$dyad-73..$dyad+73]),"\n";
	$overlapping_calls_nuc{$dyad}=$symcurv_nuc[$dyad];  #creating hash non overlapping calls
	}
}
#print " -> ",scalar(keys %overlapping_calls_nuc)," overlapping calls for 'stationary' state stored in: out_init_calls_nuc.gff\n";

# 2. "activated" state
my %overlapping_calls_dnase = ();
for ($dyad=1; $dyad<=scalar(@symcurv_dnase); $dyad++) {
	if ((defined $symcurv_dnase[$dyad]) and ($symcurv_dnase[$dyad]>0) and ($symcurv_dnase[$dyad]< 100) and ($dyad-73>0) and ($dyad+73 < length($sequence))) {
	print OUT3 $title,"\tevidence\t","act_nucleosome\t",$dyad-73+$offset,"\t",$dyad+73+$offset,"\t",$symcurv_dnase[$dyad],"\t+\t.\t",join('', @DNA[$dyad-73..$dyad+73]),"\n";
	$overlapping_calls_dnase{$dyad}=$symcurv_dnase[$dyad];  #creating hash of overlapping calls
	}
	if ((defined $symcurv_dnase[$dyad]) and ($symcurv_nuc[$dyad]>=100) and ($dyad-73>0) and ($dyad+73 < length($sequence))) {
	print OUT3 $title,"\tevidence\t","act_nucleosome\t",$dyad-73+$offset,"\t",$dyad+73+$offset,"\t100.0\t+\t.\t",join('', @DNA[$dyad-73..$dyad+73]),"\n";
	$overlapping_calls_dnase{$dyad}=$symcurv_dnase[$dyad];  #creating hash of overlapping calls
	}
}
#print " -> ",scalar(keys %overlapping_calls_dnase)," overlapping calls for 'activated' state stored in: out_init_calls_dnase.gff\n";


# 
# Calling Greedy algorithm subroutine for non overlapping positioning
# 1. For "stationary" state
my %non_overlapping_calls_nuc = GREEDYPOS(%overlapping_calls_nuc);
%non_overlapping_calls_nuc = %non_overlapping_calls;
#
foreach $dyad (sort { $a <=> $b } keys %non_overlapping_calls_nuc) {
	print OUT4 $title,"\tevidence\t","stat_nucleosome\t",$dyad-73+$offset,"\t",$dyad+73+$offset,"\t",$non_overlapping_calls_nuc{$dyad},"\t+\t.\t",join('', @DNA[$dyad-73..$dyad+73]),"\n";
	}
#print " -> ",scalar(keys %non_overlapping_calls_nuc)," non overlapping calls for 'stationary' state stored in: out_fin_calls_nuc.gff\n";

# 2. For "activated" state
my %non_overlapping_calls_dnase = GREEDYPOS(%overlapping_calls_dnase);
#
%non_overlapping_calls_dnase = %non_overlapping_calls;
foreach $dyad (sort { $a <=> $b } keys %non_overlapping_calls_dnase) {
	print OUT5 $title,"\tevidence\t","act_nucleosome\t",$dyad-73+$offset,"\t",$dyad+73+$offset,"\t",$non_overlapping_calls_dnase{$dyad},"\t+\t.\t",join('', @DNA[$dyad-73..$dyad+73]),"\n";
	}
#print " -> ",scalar(keys %non_overlapping_calls_dnase)," non overlapping calls for 'activated' state stored in: out_fin_calls_dnase.gff\n";

# 	 		}# checking sequence size
#		else {
#		$average=0;
#		}
	}#closing if 'end of fasta seq'.penultimate loop	
#closing ARGV reading. Ultimate loop
}

# closing of Curvature file
close OUT0;
# closing SymCurv file
close OUT1;
#
close OUT2;
close OUT3;
#
close OUT4;
close OUT5;

$stop = time (); #stopping timer;
print "Run complete in ",$stop-$start," seconds\n";

#End of Program here
exit ;


####SUBROUTINES#########
##########################################################################################################################
#-------------------------------
# READFASTA: Reading Sequence subroutine 
sub READFASTA{
#
undef(@DNA); #initialization
undef(@NumDNA); #initialization
# getting rid of spaces,blanks and numbering
$sequence =~ s/\n//g;
$sequence =~ tr/a-z/A-Z/;
#Caution: Ns are changed into As
$sequence =~ s/[^AGCT]/A/g;
#Case: "plus" strand -> forward 
if ($strand eq "+") {
@DNA = (split ('',$sequence));
}
#Case "minus" strand -> revcomp 
if ($strand eq "-") {
$sequence =~ tr/AGCT/TCGA/;
@DNA = reverse (split ('',$sequence));
}
for (my $i=0; $i <= scalar(@DNA); $i++) {
	if ($DNA[$i] eq "A") {$NumDNA[$i] = "0"}
	if ($DNA[$i] eq "T") {$NumDNA[$i] = "1"}
	if ($DNA[$i] eq "G") {$NumDNA[$i] = "2"}
	if ($DNA[$i] eq "C") {$NumDNA[$i] = "3"}
	}
return @DNA;
return @NumDNA;
}
#-------------------------------



#Curvature Subroutine 
sub CURVATURE {
#
#Setting structural parameters in arrays:
#Curvature parameters as published by Munteanu et al., 1998 
#Roll angles are the definitive parameter. Subscript "nuc" means nucleosomes, subscript "simple" means DNase.
#
@twist = (
[
	["0.598647428", "0.598647428", "0.598647428", "0.598647428"],
	["0.598647428", "0.598647428", "0.598647428", "0.598647428"],
	["0.598647428", "0.598647428", "0.598647428", "0.598647428"],
	["0.598647428", "0.598647428", "0.598647428", "0.598647428"],
],
[
	["0.598647428", "0.598647428", "0.598647428", "0.598647428"],
	["0.598647428", "0.598647428", "0.598647428", "0.598647428"],
	["0.598647428", "0.598647428", "0.598647428", "0.598647428"],
	["0.598647428", "0.598647428", "0.598647428", "0.598647428"],
],
[
	["0.598647428", "0.598647428", "0.598647428", "0.598647428"],
	["0.598647428", "0.598647428", "0.598647428", "0.598647428"],
	["0.598647428", "0.598647428", "0.598647428", "0.598647428"],
	["0.598647428", "0.598647428", "0.598647428", "0.598647428"],
],
[
	["0.598647428", "0.598647428", "0.598647428", "0.598647428"],
	["0.598647428", "0.598647428", "0.598647428", "0.598647428"],
	["0.598647428", "0.598647428", "0.598647428", "0.598647428"],
	["0.598647428", "0.598647428", "0.598647428", "0.598647428"],
],
);
@roll_nuc = (
[
	["0.0633", "0.3500", "4.6709", "2.64115"],
	["6.2734", "0.3500", "7.7171", "4.44325"],
	["4.8884", "3.9232", "5.0523", "6.8829"],
	["5.4903", "3.9232", "5.3055", "5.3055"],
],
[
	["4.6709", "6.2734", "5.00295", "5.0673"],
	["4.6709", "0.0633", "4.7618", "4.0633"],
	["7.7000", "5.4903", "3.05865", "6.75525"],
	["7.7000", "4.8884", "7.07195", "4.9907"],
],
[
	["4.0633", "4.44325", "5.9806", "5.51645"],
	["5.0673", "2.64115", "6.62555", "5.51645"],
	["4.9907", "5.3055", "5.89135", "9.0823"],
	["6.75525", "6.8829", "5.89135", "9.0823"],
],
[
	["4.7618", "7.7171", "6.8996", "6.62555"],
	["5.00295", "4.6709", "6.8996", "5.9806"],
	["7.07195", "5.3055", "3.869", "5.9000"],
	["3.05865", "5.0523", "3.869", "5.827"],
],
);
@roll_simple = (
[
	["0.1", "0.0", "4.2", "1.6"],
	["9.7", "0.0", "8.7", "3.6"],
	["6.5", "2.0", "4.7", "6.3"],
	["5.8", "2.0", "5.2", "5.2"],
],
[
	["7.3", "9.7", "7.8", "6.4"],
	["7.3", "0.1", "6.2", "5.1"],
	["10.0", "5.8", "0.7", "7.5"],
	["10.0", "6.5", "5.8", "6.2"],
],
[
	["5.1", "3.6", "6.6", "5.6"],
	["6.4", "1.6", "6.8", "5.6"],
	["6.2", "5.2", "5.7", "8.2"],
	["7.5", "6.3", "4.3", "8.2"],
],
[
	["6.2", "8.7", "9.6", "6.8"],
	["7.8", "4.2", "9.6", "6.6"],
	["5.8", "5.2", "3.0", "4.3"],
	["0.7", "4.7", "3.0", "5.7"],
],
);
@tilt = (
[
	["0.0", "0.0", "0.0", "0.0"],
	["0.0", "0.0", "0.0", "0.0"],
	["0.0", "0.0", "0.0", "0.0"],
	["0.0", "0.0", "0.0", "0.0"],
],
[
	["0.0", "0.0", "0.0", "0.0"],
	["0.0", "0.0", "0.0", "0.0"],
	["0.0", "0.0", "0.0", "0.0"],
	["0.0", "0.0", "0.0", "0.0"],
],
[
	["0.0", "0.0", "0.0", "0.0"],
	["0.0", "0.0", "0.0", "0.0"],
	["0.0", "0.0", "0.0", "0.0"],
	["0.0", "0.0", "0.0", "0.0"],
],
[
	["0.0", "0.0", "0.0", "0.0"],
	["0.0", "0.0", "0.0", "0.0"],
	["0.0", "0.0", "0.0", "0.0"],
	["0.0", "0.0", "0.0", "0.0"],
],
);
#Curvature Initializations
#
@xcoord_nuc =0;
@ycoord_nuc =0;
$rxsum_nuc =0;
$rysum_nuc =0;
@xave_nuc =0;
@yave_nuc =0;
@xcoord_simple =0;
@ycoord_simple =0;
$rxsum_simple =0;
$rysum_simple =0;
@xave_simple =0;
@yave_simple =0;
$twistsum = 0;
$pi = 3.14159;
undef (@Curv_nuc); #crucial for multi fasta files. Array needs to be re-initialized
undef (@Curv_simple);#crucial for multi fasta files. Array needs to be re-initialized
#default curvstep is 15
$curvstep = 15;
$curvscale = 0.33335;
#
# getting the DNA numerical array
@CurvDNA = @_;
#
# calculating curvature
for (my $j = 0; $j < scalar(@CurvDNA)-2; $j++){
	$twistsum += $twist[$CurvDNA[$j]][$CurvDNA[$j+1]][$CurvDNA[$j+2]];
	my $dx_nuc = $roll_nuc[$CurvDNA[$j]][$CurvDNA[$j+1]][$CurvDNA[$j+2]]*sin($twistsum)+ $tilt[$CurvDNA[$j]][$CurvDNA[$j+1]][$CurvDNA[$j+2]]*sin($twistsum-$pi/2.0);
	my $dx_simple = $roll_simple[$CurvDNA[$j]][$CurvDNA[$j+1]][$CurvDNA[$j+2]]*sin($twistsum)+ $tilt[$CurvDNA[$j]][$CurvDNA[$j+1]][$CurvDNA[$j+2]]*sin($twistsum-$pi/2.0);
	my $dy_nuc = $roll_nuc[$CurvDNA[$j]][$CurvDNA[$j+1]][$CurvDNA[$j+2]]*cos($twistsum)+ $tilt[$CurvDNA[$j]][$CurvDNA[$j+1]][$CurvDNA[$j+2]]*cos($twistsum-$pi/2.0);
	my $dy_simple = $roll_simple[$CurvDNA[$j]][$CurvDNA[$j+1]][$CurvDNA[$j+2]]*cos($twistsum)+ $tilt[$CurvDNA[$j]][$CurvDNA[$j+1]][$CurvDNA[$j+2]]*cos($twistsum-$pi/2.0);
	$xcoord_nuc[$j+1]=$xcoord_nuc[$j]+$dx_nuc;
	$ycoord_nuc[$j+1]=$ycoord_nuc[$j]+$dy_nuc;
	$xcoord_simple[$j+1]=$xcoord_simple[$j]+$dx_simple;
	$ycoord_simple[$j+1]=$ycoord_simple[$j]+$dy_simple;
	}
#additional parameters for curvature calculation
#default values $stepone=6, $steptwo=4
my $stepone=6.0;
my $steptwo=4.0;
for (my $k = $stepone; $k < scalar(@CurvDNA)-$stepone; $k++) {
		$rxsum_nuc = 0;
		$rysum_nuc = 0;
		$rxsum_simple = 0;
		$rysum_simple = 0;
	for ($ik = -$steptwo; $ik <= $steptwo; $ik++) {
		$rxsum_nuc += $xcoord_nuc[$k+$ik];
		$rysum_nuc += $ycoord_nuc[$k+$ik];
		$rxsum_simple += $xcoord_simple[$k+$ik];
		$rysum_simple += $ycoord_simple[$k+$ik];
		}
		$rxsum_nuc += $xcoord_nuc[$k+($stepone-1)]/($steptwo/2);
		$rysum_nuc += $ycoord_nuc[$k+($stepone-1)]/($steptwo/2);
		$rxsum_nuc += $xcoord_nuc[$k-($stepone-1)]/($steptwo/2);
		$rysum_nuc += $ycoord_nuc[$k-($stepone-1)]/($steptwo/2);
		$xave_nuc[$k] = $rxsum_nuc/(($stepone-1)*2);
		$yave_nuc[$k] = $rysum_nuc/(($stepone-1)*2);
		$rxsum_simple += $xcoord_simple[$k+($stepone-1)]/($steptwo/2);
		$rysum_simple += $ycoord_simple[$k+($stepone-1)]/($steptwo/2);
		$rxsum_simple += $xcoord_simple[$k-($stepone-1)]/($steptwo/2);
		$rysum_simple += $ycoord_simple[$k-($stepone-1)]/($steptwo/2);
		$xave_simple[$k] = $rxsum_simple/(($stepone-1)*2);
		$yave_simple[$k] = $rysum_simple/(($stepone-1)*2);
	}
for (my $l = $curvstep+$stepone; $l < scalar(@CurvDNA)-$curvstep-$stepone; $l++){
	$Curv_nuc[$l] = (sqrt( ($xave_nuc[$l+$curvstep]-$xave_nuc[$l-$curvstep])**2 + ($yave_nuc[$l+$curvstep]-$yave_nuc[$l-$curvstep])**2 ));
	$Curv_simple[$l] = (sqrt( ($xave_simple[$l+$curvstep]-$xave_simple[$l-$curvstep])**2 + ($yave_simple[$l+$curvstep]-$yave_simple[$l-$curvstep])**2 ));
	$Curv_nuc[$l] *= $curvscale;
	$Curv_simple[$l] *= $curvscale;
	}
#
return @Curv_nuc;
return @Curv_simple;
}
##########################################################################################################################


##########################################################################################################################
# SymCurv subroutine
sub SYMCURV {
#
@Curv_nuc = @_;
$average=0;
undef (@symcurv); #crucial for multi fasta files. Array needs to be re-initialized
#default window size is 51;
my $win = 101;
my $step = 1;
#Part 1: breaking up in core pieces of 146 bps
	for (my $dyad = $win ; $dyad < scalar(@Curv_nuc) - $win ; $dyad += $step) {
		my $weight = 0;
		my $sum = 0;
#Part 2: symmetry calculations
		for (my $j = $dyad , my $k = $dyad ; $j < $dyad + int($win/2) + 1, $k > $dyad - int($win/2) - 1 ; $j += $step, $k -= $step) {
		#for (my $j = $dyad  ; $j < $dyad + int($win/2) + 1, ; $j += $step) {
			#Calculating symmetry of curvature
			#Calculating mirror coordinate in respect to $dyad (if x is dyad mirror is 2*x-j)
			$score = abs($Curv_nuc[$j] - $Curv_nuc[$k]);
			$sum += $score;
			}
				if ( ($Curv_nuc[$dyad] < $Curv_nuc[$dyad-1]) and ($Curv_nuc[$dyad] < $Curv_nuc[$dyad+1]) and ( (($Curv_nuc[$dyad-1] - $Curv_nuc[$dyad]) + ($Curv_nuc[$dyad+1] - $Curv_nuc[$dyad])) >= 0.01 )) {
				$weight =  1/(($Curv_nuc[$dyad-1] - $Curv_nuc[$dyad]) + ($Curv_nuc[$dyad+1] - $Curv_nuc[$dyad]))
				}
			# Minimum set to 0 instead of 1 when treating non minima
				else {
				$weight = 0;
				}
					#final calculation of symcurv values
					if ($sum != 0) {
					$symcurv[$dyad] = (1/$sum) * $weight;
					$average += $symcurv[$dyad];
					}
					else {
					print "Symmetry component is 0 for position ", $dyad+$offset,"\n";
					$symcurv[$dyad] = 100;
					$average += $symcurv[$dyad];
					}
			}
$average/=scalar(@Curv_nuc);		
return @symcurv;
return $average;
}
##########################################################################################################################


##########################################################################################################################
sub GREEDYPOS {
#initializations
my $spacer = 30 ; # minimum size of the space between two consecutive non-overlapping segments
my $size = 147 ; #size of desired non-overlapping segments. if variable it can be generated by a function
my $num = 0 ;
my $no = 0 ;
my $pos = 0 ; #intermediate position scalar
my $score = 0; #intermediate score scalar;
my @position = 0; #final positions array
my @score_nuc = 0 ; #final score array
%overlapping_calls=@_;
undef (%non_overlapping_calls); #crucial for multi fasta files. Hash needs to be re-initialized
#sorting hash in descending order by value
foreach $dyad (sort { $overlapping_calls{$b} <=> $overlapping_calls{$a} } keys %overlapping_calls) {
		$num++;
		$pos = $dyad;
		$score= $overlapping_calls{$dyad};
		$index = 0;
		for (my $i = 0; $i <= $no; $i++) {
		# testing for the negative case, i.e. segment overlapping one of the previous segments
			if ( ($pos >= $position[$i]-$spacer-$size) and ($pos <= $position[$i]+$spacer+$size)) { 
			$index = 1;
			}
		}
		if ($index == 0) {
		$no++; 
		$position[$no] = $pos; 
		$score_nuc[$no] = $score;
		}
     	}
     	# writing output in a hash to return as product of the subroutine
	for (my $j = 1; $j <= $no; $j++) {
	$non_overlapping_calls{$position[$j]}=$score_nuc[$j];
	}
return %non_overlapping_calls;
}
##########################################################################################################################
