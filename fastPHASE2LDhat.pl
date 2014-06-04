#! /usr/bin/perl

#program converts fastPHASE output to LDhat sites input
#locs input file should be generated with 'vcf2fastPHASE.pl' program
#August 13th, 2012
#Modified September 4, 2012

$fP_output = $ARGV[0];
$fP_input = $ARGV[1];
$positions = $ARGV[2];

unless ($#ARGV==2) {
    print STDERR "Please provide name of fastPHASE output filename, fastPHASE input filename, and the positions filename on command line\n\n";
    die;
} #end unless

print STDERR "Reading in positions file...";
open(POS, $positions);
@positions_start = ();
@positions_end = ();

while(<POS>) {
    chomp;
    @input = split(/\s+/, $_);
    push(@positions_start, $input[0]);
    push(@positions_end, $input[1]);
#    print STDERR "\n$input[0]\t$input[1]\n"; 
} #end while

print STDERR "done.\nReading in fastPHASE output file...";
open(INPUT, $fP_output);
$start_input = 0;
$start_new_indv = 0;
$start_haplotype1 = 0;
$start_haplotype2 = 0;

%haplotypes = ();
@hap_names = ();

while(<INPUT>) {
    chomp;
    if ($_=~/BEGIN GENOTYPES/) {
	$start_input = 1;
	$start_new_indv = 1;
	next;
    } elsif ($start_input==0) {
	next;
    } elsif ($_=~/END GENOTYPES/) {
	last;
    } #end elsif

    if ($start_new_indv==1 && $start_haplotype1==0) {
	@input_array = split(/\s+/, $_);
	$hap1_name = ">" . $input_array[1] . "-0";
	$hap2_name = ">" . $input_array[1] . "-1";
#	print STDERR "$hap1_name\n";
	push(@hap_names, $hap1_name);
	$start_haplotype1 = 1;
	$start_new_indv = 0;
    } elsif ($start_haplotype1==1 && $start_haplotype2==0) {
	@haplotype = split(/\s+/, $_);
	for ($i=0; $i<=$#haplotype; $i++) {
	    push @{$haplotypes{$hap1_name}}, $haplotype[$i];
	} #end for
	push(@hap_names, $hap2_name);
#	print STDERR "$hap2_name\n";
	$start_haplotype2 = 1;
	$start_haplotype1 = 0;
    } elsif ($start_haplotype1==0 && $start_haplotype2==1) {
	@haplotype = split(/\s+/, $_);
	for ($i=0; $i<=$#haplotype; $i++) {
	    push @{$haplotypes{$hap2_name}}, $haplotype[$i];
	} #end for
	$start_haplotype2 = 0;
	$start_new_indv = 1;
    } #end else

} #end while

print STDERR "done.\nWriting LDhat output files...";

@fP_pieces1 = split(/\./, $fP_input);
#print STDERR "fastPHASE input: $fP_input; chr: $fP_pieces1[0]; block: $fP_pieces1[2]\n";
#@fP_pieces2 = split(/\./, $fP_pieces1[$#fP_pieces1]);

$positions_file = "$fP_pieces1[0].synteny_block.$fP_pieces1[2].positions.out";

open(POS2, $positions_file);
$header = (<POS2>);
@sites=();

while (<POS2>) {
    chomp;
    push(@sites, $_);
} #end while                    

open(fastPHASE, $fP_input);
$num_indv = (<fastPHASE>);
chomp($num_indv);
$num_hap = $num_indv*2;
$num_sites = (<fastPHASE>);
chomp($num_sites);
#$sites_input = (<fastPHASE>);
#chomp($sites_input);
#@sites = split(/\s+/, $sites_input);
#shift(@sites);

#print STDERR "Length of positions array: $#positions_start\n";

for ($i=0; $i<=$#positions_start; $i++) {

    $count = $i+1;
    $locs = "LDhat-inputs\/$fP_pieces1[0]\/$fP_pieces1[0]\.$positions_start[$i]\.$positions_end[$i]\.ldhat\.locs";
    open(LOCS, ">$locs");
	
    for ($l=0; $l<=$#sites; $l++) {
	if ($sites[$l]==$positions_start[$i]) {
	    $start_k = $l;
	} elsif ($sites[$l]==$positions_end[$i]) {
	    $end_k = $l;
	    last;
	} #end elsif
    } #end for

    $num_sites_i_block = $end_k - $start_k + 1;
    $end_block = $sites[$end_k]/1000;
    print LOCS "$num_sites_i_block $end_block L\n";    

    $sites = "LDhat-inputs\/$fP_pieces1[0]\/$fP_pieces1[0]\.$positions_start[$i]\.$positions_end[$i]\.ldhat\.sites";
    open(SITES, ">$sites");
    print SITES "$num_hap $num_sites_i_block 1";

    for ($j=0; $j<=$#hap_names; $j++) {
	$current_hap_name = $hap_names[$j];
	@current_hap = @{$haplotypes{$current_hap_name}};
	print SITES "\n$current_hap_name\n";
	for ($k=$start_k; $k<=$end_k; $k++) {
	    print SITES "$current_hap[$k]";
	    if ($j==0) {
		$kb_sites = $sites[$k]/1000;
		printf LOCS ("%0.3f\n", $kb_sites);
	    } #end if
	} #end for
    } #end for

    $list = ">>$fP_pieces1[0]-input.BED";
#    print STDERR "cat \"$fP_pieces1[0]\t$sites[$start_k]\t$sites[$end_k]\" $list\n";
    system("echo \"$fP_pieces1[0]\t$sites[$start_k]\t$sites[$end_k]\" $list");
} #end for

print STDERR "done.\n";
