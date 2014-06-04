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
    $ref = $input_line[3];
    $alt = $input_line[4];

    for ($i=9; $i<=$loop_size; $i++) {
	$o = $i - 9;
	$hap1 = $names[$o];
	@genotype = split(":", $input_line[$i]);
#	print STDERR "Hap name: $hap1; Original genotype: $input_line[$i]; Next step: $genotype[0];";
	$genotype[0] =~ s/0/$ref/g;
	$genotype[0] =~ s/1/$alt/g;
	$genotype[0] =~ s/\./\?/g;
#	print STDERR " Final: $genotype[0]\n";
	push @{$genotypes{$hap1}}, $genotype[0];
    } #end for
} #end while

print STDERR "done.\nNow printing output...";

open(OUTPUT, ">$output");
$hap_size = $sample_size;
print OUTPUT "$hap_size";

for ($c=0; $c<=$#names; $c++) {
    print OUTPUT "\n\# $names[$c]\n";

    $hap1 = $names[$c];
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

#    print OUTPUT "\n$hap1_line";

} #end for

print OUTPUT "\n";

print STDERR "done.\n";
