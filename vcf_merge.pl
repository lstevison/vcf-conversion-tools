#! /usr/bin/perl

#program combined two vcf files with different sets of inidividuals into one file
#July 15, 2013
#Last Modified July 15, 2013

$VCF1 = $ARGV[0];
$VCF2 = $ARGV[1];
$VCF_out = $ARGV[2]; 

unless ($#ARGV==2) {
    print STDERR "Please provide input names of both VCF files and output filename on command line.\n\n";
    die;
} #end unless

print STDERR "Now reading through two VCF input files and writing merged output...\n";

open(OUTPUT, ">$VCF_out");
open(VCF1, $VCF1);
open(VCF2, $VCF2);

while (<VCF1>) {
    chomp;
    if ($_=~/\#\#/) {
	print OUTPUT "$_\n";
        next;
    } elsif ($_=~/\#CHROM/) {

	 print OUTPUT "$_";

	while (<VCF2>) {
	    chomp;
	    if ($_=~/\#CHROM/) {
		@vcf2_chrom_line = split(/\s+/, $_);
		last;
	    } #end if
	} #end while

	for ($d=9; $d<=$#vcf2_chrom_line; $d++) {
	    print OUTPUT "\t$vcf2_chrom_line[$d]";
	} #end for

        next;
    } elsif ($_!~/chr/) {
	next;
    } #end elsif

    @input_line = split(/\s+/, $_);
    $vcf1_position = $input_line[1];

    $vcf2_line = (<VCF2>);
    @vcf2_lines = split(/\s+/, $vcf2_line);
    $vcf2_position = $vcf2_lines[1];

    #compare position to @original_positions array from fastPHASE input file!
    if ($vcf2_position!=$vcf1_position) {
	print STDERR "Error: position in VCF file 1 ($vcf1_position) does not match position in second VCF file ($vcf2_position)\n";
	next;
    } #end if

    print OUTPUT "\n$_";

    for ($e=9; $e<=$#vcf2_lines; $e++) {
	print OUTPUT "\t$vcf2_lines[$e]";
    } #end for


} #end while                           


print OUTPUT "\n";
print STDERR "done.\n";
