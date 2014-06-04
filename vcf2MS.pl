#! /usr/bin/perl

#program converts vcf file to MS input file format
#January 17th, 2013

$vcf = $ARGV[0];
$output = $ARGV[1];
$sample_size = $ARGV[2];

unless ($#ARGV==2) {
    print STDERR "Please provide name of input vcf file, filename for MS formatted output, and sample size on command line\n\n";
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

#    $ref = $input_line[3];
#    $alt = $input_line[4];

    for ($i=9; $i<=$loop_size; $i++) {
	$o = $i - 9;
	$hap1 = $names[$o] . "_1";
	$hap2 = $names[$o] . "_2";
	@genotype = split(":", $input_line[$i]);
#	print STDERR "Hap1: $hap1; Hap2: $hap2; $Before: $genotype[0]\t";
#	$genotype[0] =~ s/0/$ref/g;
#	$genotype[0] =~ s/1/$alt/g;
#	$genotype[0] =~ s/\./\?/g;
#	print STDERR "After: $genotype[0]\n";
	@haplotypes = split(/[\|\/]/, $genotype[0]);
	push @{$genotypes{$hap1}}, $haplotypes[0];
	push @{$genotypes{$hap2}}, $haplotypes[1];
    } #end for
} #end while

print STDERR "done.\nNow printing output...";

open(OUTPUT, ">$output");
$number_loci = $#positions + 1;

print OUTPUT "$number_loci\n";

for ($b=0; $b<=$#positions; $b++) {
	print OUTPUT "$positions[$b] ";
} #end for

for ($c=0; $c<=$#names; $c++) {
    print OUTPUT "\n";

    $hap1 = $names[$c] . "_1";
    @hap1_geno = @{$genotypes{$hap1}};

    for ($d=0; $d<=$#hap1_geno; $d++) {
	    print OUTPUT "$hap1_geno[$d]";
    } #end for
    print OUTPUT "\n";

    $hap2 = $names[$c] . "_2";
    @hap2_geno = @{$genotypes{$hap2}};

    for ($e=0; $e<=$#hap2_geno; $e++) {
	    print OUTPUT "$hap2_geno[$e]";
    } #end for
} #end for

print OUTPUT "\n";

print STDERR "done.\n";
