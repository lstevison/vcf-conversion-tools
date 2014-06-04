#! /usr/bin/perl

#program converts vcf file to fastPHASE input format
#September 13th, 2012

$vcf = $ARGV[0];
$output = $ARGV[1];
$output2 = $ARGV[2];
$sample_size = $ARGV[3];

unless ($#ARGV==3) {
    print STDERR "Please provide name of input vcf file, filename for fastPHASE formatted output, filename for positions file, and sample size on command line\n\n";
    die;
} #end unless

open(VCF, $vcf);

@positions = ();
@names = ();
$loop_size = $sample_size + 8;
%genotypes = ();

print STDERR "Reading in VCF file...";

while(<VCF>) {
    chomp;
    if ($_=~/\#\#/) {
	next;
    } elsif ($_=~/\#/) {
	@input_line = split(/\s+/, $_);
	for ($a=9; $a<=$loop_size; $a++) {
	    push(@names, $input_line[$a]);
	} #end for
	next;
    } #end elsif
    
    @input_line = split(/\s+/, $_);
    push(@positions, $input_line[1]);

    $ref = $input_line[3];
    $alt = $input_line[4];

    for ($i=9; $i<=$loop_size; $i++) {
	$o = $i - 9;
	$hap1 = $names[$o] . "_1";
	$hap2 = $names[$o] . "_2";
	@genotype = split(":", $input_line[$i]);
#	print STDERR "Hap1: $hap1; Hap2: $hap2; $Before: $genotype[0]\t";
	$genotype[0] =~ s/0/$ref/g;
	$genotype[0] =~ s/1/$alt/g;
	$genotype[0] =~ s/\./\?/g;
#	print STDERR "After: $genotype[0]\n";
	@haplotypes = split(/[\|\/]/, $genotype[0]);
	push @{$genotypes{$hap1}}, $haplotypes[0];
	push @{$genotypes{$hap2}}, $haplotypes[1];
    } #end for
} #end while

print STDERR "done.\nNow printing output...";

open(OUTPUT, ">$output");
open(OUTPUT2, ">$output2");

$number_loci = $#positions + 1;

print OUTPUT "$sample_size\n$number_loci\nP ";
print OUTPUT2 "CHR\tPOS\n";
$positions_line = "P ";

for ($b=0; $b<=$#positions; $b++) {
    $positions_line .= "$positions[$b] ";
    $cnt = length($positions_line);

    if ($cnt<500000) {
	print OUTPUT "$positions[$b] ";
    } elsif ($cnt>=500000) {
	print OUTPUT "\n$positions[$b] ";
	$positions_line = "$positions[$b] ";
    } #end elsif
    print OUTPUT2 "$positions[$b]\n";
} #end for

for ($c=0; $c<=$#names; $c++) {
    print OUTPUT "\n\# $names[$c]\n";

    $hap1 = $names[$c] . "_1";
    @hap1_geno = @{$genotypes{$hap1}};
    $hap1_line = "";

    for ($d=0; $d<=$#hap1_geno; $d++) {
	$hap1_line .= $hap1_geno[$d];
	$cnt2 = length($hap1_line);
	if ($cnt2<500000) {
	    print OUTPUT "$hap1_geno[$d]";
	} elsif ($cnt2>=500000) {
	    print OUTPUT "\n$hap1_geno[$d]";
	    $hap1_line = "";
	} #end elsif
    } #end for
    print OUTPUT "\n";

    $hap2 = $names[$c] . "_2";
    @hap2_geno = @{$genotypes{$hap2}};
    $hap2_line = "";

    for ($e=0; $e<=$#hap2_geno; $e++) {
	$hap2_line .= $hap2_geno[$e];
	$cnt3 = length($hap2_line);
	if ($cnt3<500000) {
	    print OUTPUT "$hap2_geno[$e]";
	} elsif ($cnt3>=500000) {
	    print OUTPUT "\n$hap2_geno[$e]";
	    $hap2_line = "";
	} #end elsif
    } #end for
} #end for

print OUTPUT "\n";

print STDERR "done.\n";
