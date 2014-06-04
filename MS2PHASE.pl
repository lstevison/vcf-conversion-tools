#! /usr/bin/perl

#Program reads in MS formatted input (*.hud) and converts to PHASE input format

$input=$ARGV[0];
$output=$ARGV[1];
$num_ind=$ARGV[2];
$num_haps=($num_ind*2) + 1;

unless ($#ARGV==2) {
    print STDERR "Please provide input and output filenames and number of individuals on command line.\n\n";
    die;
} #end unless 

open(INPUT, $input);
open(OUTPUT, ">$output");

print OUTPUT "$num_ind\n";
$loop_counter=0;

while (<INPUT>) {
    chomp;
    if ($loop_counter >$num_haps) {
	last;
    } #end if

    if ($loop_counter==0) {
	$num_loci=$_;
	print OUTPUT "$num_loci\n";    #prints first line of input into output (number of loci)
	$loop_counter++;
	next;
    } elsif ($loop_counter==1) {       #prints next two lines of output with positions and specifies all loci to be SNPs
	@positions=split(/\s+/,$_);
	print OUTPUT "P";
	for ($i=0; $i<@positions; $i++) {
	    $locus=$positions[$i]*1000000;
	    print OUTPUT " $locus";
	} #end for

	print OUTPUT "\n";
	for ($s=0; $s<@positions; $s++) {
	    print OUTPUT "S";
	} #end for
	print OUTPUT "\n";
	$loop_counter++;
	next;
    } #end elsif

    if ($loop_counter % 2 == 0) { #loop counter is even (need to specify haplotype names)
	print OUTPUT "#$loop_counter\n$_\n";
    } else {
	print OUTPUT "$_\n";
    } #end else
    
    $loop_counter++;

} #end while
