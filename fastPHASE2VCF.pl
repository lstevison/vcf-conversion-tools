#! /usr/bin/perl

#program converts fastPHASE output to VCF file
#program needs original VCF file and input file for fastPHASE
#January 4th, 2013
#Last Modified Jan 8th, 2013

$fP_input = $ARGV[0];
$fP_output = $ARGV[1];
$VCF = $ARGV[2];
$VCF_out=$ARGV[3];
$block=$ARGV[4];

@filename = split(/\./, $fP_input);
@chr = split(/\//, $filename[0]);
$error_report = "error_report.$chr[1].txt";

unless ($#ARGV==4) {
    print STDERR "Please provide name of fastPHASE input and output filenames, the original VCF filename and an output VCF filename, and the synteny block number on command line\n\n";
    die;
} #end unless

print STDERR "Reading in fastPHASE input file...";
open(FPINPUT, $fP_input);

@original_positions = ();
$positions_line = 0;

while (<FPINPUT>) {
    chomp;
    if ($_=~/P\s+/) {
	$positions_line = 1;
    } elsif ($_=~/\#/) {
	last;
    } #end elsif

    if ($positions_line==1) {
	@line_array = split(/\s+/, $_);
	push(@original_positions, @line_array);
    } #end if

} #end while

shift @original_positions;

print STDERR "positions: $original_positions[0]-$original_positions[$#original_positions]...done.\nReading in fastPHASE output file...";
open(FPOUTPUT, $fP_output);
$start_input = 0;
$start_new_indv = 0;
$start_haplotype1 = 0;
$start_haplotype2 = 0;

%haplotypes = ();
@hap_names = ();
@ind_names = ();

while(<FPOUTPUT>) {
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
	$hap1_name = $input_array[1] . "_1";
	$hap2_name = $input_array[1] . "_2";
#	print STDERR "$hap1_name\n";
	push(@ind_names, $input_array[1]);
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

#compare position to @original_positions array from fastPHASE input file!
@output_size = @{$haplotypes{$hap2_name}};
if ($#output_size!=$#original_positions) {
    print STDERR "Error: Original number of sites does not match output number of sites in fastPHASE files!\n";
    die;
} #end if

$loop_size=$#ind_names + 9;
$num_indv = $#ind_names + 1;
$while_loop_counter = 0;
@names= ();
print STDERR "done.\nReading in original VCF file with $num_indv individuals and writing output VCF file...";
open(OUTPUT, ">$VCF_out");
open(VCF, $VCF);

while (<VCF>) {
    chomp;
    if ($_=~/\#\#/) {
	print OUTPUT "$_\n";
        next;
    } elsif ($_=~/\#/) {
        @input_line = split(/\s+/, $_);
        for ($a=9; $a<=$loop_size; $a++) {
	    $b = $a - 9;
	    push(@names, $input_line[$a]);
            if ($input_line[$a]!~/$ind_names[$b]/) {
		print STDERR "error: names in fastPHASE file and VCF don't match\n\n";
		die;
	    } #end if
        } #end for                                                                                                                                                                                              

	print OUTPUT "$_";
        next;
    } #end elsif                                                                                                                                                                                                 
    @input_line = split(/\s+/, $_);
    $position = $input_line[1];

    #compare position to @original_positions array from fastPHASE input file!
    if ($original_positions[$while_loop_counter]!=$position) {
	print STDERR "Error: position in VCF file ($position) does not match position in fastPHASE input file ($original_positions[$while_loop_counter])!\n";
	next;
    } #end if

    $ref = $input_line[3];
    $alt = $input_line[4];
    print OUTPUT "\n$input_line[0]\t$input_line[1]\t$input_line[2]\t$input_line[3]\t$input_line[4]\t$input_line[5]\t$input_line[6]\t$input_line[7];synteny_block=$block\t$input_line[8]";

    for ($i=9; $i<=$loop_size; $i++) {
        $o = $i - 9;
        $hap1 = $names[$o] . "_1";
        $hap2 = $names[$o] . "_2";
        @genotype = split(":", $input_line[$i]);
#       print STDERR "Hap1: $hap1; Hap2: $hap2; $Before: $genotype[0] ";                                                                                                                                       
	$allele1 = @{$haplotypes{$hap1}}[$while_loop_counter];
        $allele1 =~ s/$ref/0/g;
        $allele1 =~ s/$alt/1/g;

	$allele2 = @{$haplotypes{$hap2}}[$while_loop_counter]; 
        $allele2 =~ s/$ref/0/g;
        $allele2 =~ s/$alt/1/g;

	#genotype can either be homo-ref, homo-alt, het, or missing data
	@alleles = split(/\||\//, $genotype[0]);

	if ($alleles[0]=~/0/ && $alleles[1]=~/0/) { #original vcf is homozygous reference allele
	    if ($allele1!=0 || $allele2!=0) { #after phasing is not homo ref
		system("echo \"Potential genotyping error found at position $input_line[1] for $names[$o]; after phasing is not homo ref\" >>$error_report");
		system("echo \"Before phasing: $genotype[0]; allele 1: $alleles[0]; allele 2: $alleles[1]; After phasing: Allele 1: $allele1; Allele 2: $allele2\" >>$error_report");
	    } #end if
	} elsif ($alleles[0]=~/1/ && $alleles[1]=~/1/) { #original vcf is homozygous alternate allele
	    if ($allele1!=1 || $allele2!=1) { #after phasing not homo alt
		system("echo \"Potential genotyping error found at position $input_line[1] for $names[$o]; after phasing is not homo alt\" >>$error_report");
		system("echo \"Before phasing: $genotype[0]; allele 1: $alleles[0]; allele 2: $alleles[1]; After phasing: Allele 1: $allele1; Allele 2: $allele2\" >>$error_report");
	    } #end if
	} elsif (($alleles[0]=~/0/ && $alleles[1]=~/1/) || ($alleles[0]=~/1/ && $alleles[1]=~/0/)) { #original vcf is heterozygous
	    if (($allele1==1 && $allele2==1) || ($allele1==0 && $allele2==0)) { #after phasing is homozygous
		system("echo \"Potential genotyping error found at position $input_line[1] for $names[$o]; after phasing is homozygous\" >>$error_report");
		system("echo \"Before phasing: $genotype[0]; allele 1: $alleles[0]; allele 2: $alleles[1]; After phasing: Allele 1: $allele1; Allele 2: $allele2\" >>$error_report");
	    } #end if
	} #end elsif

	shift @genotype;
	$end_of_line = join (":", @genotype);

	print OUTPUT "\t$allele1|$allele2:$end_of_line";
    } #end for       

    $while_loop_counter++; 
} #end while                           


print OUTPUT "\n";
print STDERR "done.\n";
