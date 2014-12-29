#! /usr/bin/perl

#program converts cleaned up and combined PHASE output to VCF file
#program needs original VCF file and input file for PHASE
#January 4th, 2013
#Last Modified July 28, 2014

$P_output = $ARGV[0];
$VCF = $ARGV[1];
$VCF_out = $VCF; 
$VCF_out=~s/\.recode\.vcf/\.re-phased\.vcf/;
$block=$ARGV[2];

unless (-e $P_output) {
    print STDERR "Error: Input file does not exist: $P_output\n";
    die;
} #end if

unless (-e $VCF) {
    print STDERR "Error: Input file does not exist: $VCF\n";
    die;
} #end if

@filename = split(/(\.|\/)/, $P_output);
#@chr = split(/\./, $filename[0]);
$error_report = "error_report.$filename[2].txt";
print STDERR "Error report filename: $error_report\n";
#open(ERROR, ">>error_report.$chr[0].txt");
unless ($#ARGV==2) {
    print STDERR "Please provide name of cleaned-up and combined PHASE output file, the original VCF filename, and the synteny block number on command line\n\n";
    die;
} #end unless

open(POUTPUT, $P_output);
%haplotypes = ();
@hap_names = ();
@ind_names = ();

$num_ind=(<POUTPUT>);
chomp $num_ind;
$num_sites=(<POUTPUT>);
chomp $num_sites;
$positions=(<POUTPUT>);
chomp $positions;

$positions=~s/P //;
@original_positions=split(/\s+/,$positions);

print STDERR "Num indiv: $num_ind; Num sites: $num_sites; Positions: $original_positions[0]-$original_positions[$#original_positions]\n";

while(<POUTPUT>) {
    chomp;

    if ($_=~/-1/) {
	$start_haplotype1=1;
	$start_haplotype2=0;
	$ind_name=$_;
	$ind_name=~s/#//;
	$hap1_name=$ind_name;
	$ind_name=~s/-1//;
	push(@ind_names, $ind_name);
	push(@hap_names, $hap1_name);

	next;
    } elsif ($_=~/-2/) {
	$start_haplotype1=0;
	$start_haplotype2=1;
	$ind_name=$_;
	$ind_name=~s/#//;
	$hap2_name=$ind_name;
	$ind_name=~s/-2//;
	push(@hap_names, $hap2_name);
	next;
    } #end if

    if ($start_haplotype1==1 && $start_haplotype2==0) {
	@haplotype = split("", $_);
	for ($i=0; $i<=$#haplotype; $i++) {
	    push @{$haplotypes{$hap1_name}}, $haplotype[$i];
	} #end for
    } elsif ($start_haplotype1==0 && $start_haplotype2==1) {
	@haplotype = split("", $_);
	for ($i=0; $i<=$#haplotype; $i++) {
	    push @{$haplotypes{$hap2_name}}, $haplotype[$i];
	} #end for
    } #end else

} #end while

#print STDERR "@ind_names\n@hap_names\n";

#compare position to @original_positions array from fastPHASE input file!
@output_size = @{$haplotypes{$hap2_name}};
if ($#output_size!=$#original_positions) {
    print STDERR "Error: Original number of sites: $#original_positions does not match output number of sites in PHASE file: $#output_size!\n";
    die;
} #end if

$loop_size=$#ind_names + 9;
$num_indv = $#ind_names + 1;
$error=0;
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
    } elsif ($_=~/\#CHROM/) {
        @input_line = split(/\s+/, $_);
        for ($a=9; $a<=$loop_size; $a++) {
	    $b = $a - 9;
	    push(@names, $input_line[$a]);
            if ($input_line[$a]!~/$ind_names[$b]/) {
		print STDERR "error: names in PHASE file and VCF don't match: $input_line[$a], $ind_names[$b]\n\n";
		die;
	    } #end if
        } #end for                                                                                                                                                                                              

	print OUTPUT "$_";
        next;
    } elsif ($_!~/chr/) {
	next;
    } #end elsif

    @input_line = split(/\s+/, $_);
    $position = $input_line[1];

    #compare position to @original_positions array from PHASE input file!
    if ($original_positions[$while_loop_counter]!=$position) {
	print STDERR "Error: position in VCF file ($position) does not match position in PHASE input file ($original_positions[$while_loop_counter])!\n";
	next;
    } #end if

    $ref = $input_line[3];
    $alt = $input_line[4];
    print OUTPUT "\n$input_line[0]\t$input_line[1]\t$input_line[2]\t$input_line[3]\t$input_line[4]\t$input_line[5]\t$input_line[6]\t$input_line[7]\t$input_line[8]";

    for ($i=9; $i<=$loop_size; $i++) {
        $o = $i - 9;
        $hap1 = $names[$o] . "-1";
        $hap2 = $names[$o] . "-2";
        @genotype = split(":", $input_line[$i]);
 #      print STDERR "Hap1: $hap1; Hap2: $hap2; Before: $genotype[0] ";                                                                                                                                       
	$allele1 = @{$haplotypes{$hap1}}[$while_loop_counter];
        $allele1 =~ s/$ref/0/g;
        $allele1 =~ s/$alt/1/g;

	$allele2 = @{$haplotypes{$hap2}}[$while_loop_counter]; 
        $allele2 =~ s/$ref/0/g;
        $allele2 =~ s/$alt/1/g;

	#genotype can either be homo-ref, homo-alt, het, or missing data
#	@alleles = split(/(\/|\|)/, $genotype[0]);
	@alleles = split(/\|/, $genotype[0]);

	if ($alleles[0]=~/0/ && $alleles[1]=~/0/) { #original vcf is homozygous reference allele
	    if ($allele1!=0 || $allele2!=0) { #after phasing is not homo ref
		system("echo \"Potential genotyping error found at Block: $block position $input_line[1] for $names[$o]; after phasing is not homo ref\" >>$error_report");
		system("echo \"Before phasing: $genotype[0]; allele 1: $alleles[0]; allele 2: $alleles[1]; After phasing: Allele 1: $allele1; Allele 2: $allele2\" >>$error_report");
		$error=1;
	    } #end if
	} elsif ($alleles[0]=~/1/ && $alleles[1]=~/1/) { #original vcf is homozygous alternate allele
	    if ($allele1!=1 || $allele2!=1) { #after phasing not homo alt
		system("echo \"Potential genotyping error found at Block: $block position $input_line[1] for $names[$o]; after phasing is not homo alt\" >>$error_report");
		system("echo \"Before phasing: $genotype[0]; allele 1: $alleles[0]; allele 2: $alleles[1]; After phasing: Allele 1: $allele1; Allele 2: $allele2\" >>$error_report");
		$error=1;
	    } #end if
	} elsif (($alleles[0]=~/0/ && $alleles[1]=~/1/) || ($alleles[0]=~/1/ && $alleles[1]=~/0/)) { #original vcf is heterozygous
	    if (($allele1==1 && $allele2==1) || ($allele1==0 && $allele2==0)) { #after phasing is homozygous
		system("echo \"Potential genotyping error found at Block: $block position $input_line[1] for $names[$o]; after phasing is homozygous\" >>$error_report");
		system("echo \"Before phasing: $genotype[0]; allele 1: $alleles[0]; allele 2: $alleles[1]; After phasing: Allele 1: $allele1; Allele 2: $allele2\" >>$error_report");
		$error=1;
	    } #end if
	} #end elsif

	shift @genotype;
	$end_of_line = join (":", @genotype);
	if ($error==0) {
	    print OUTPUT "\t$allele1|$allele2:$end_of_line";
	} elsif ($error==1) {
#	    print STDERR "Error at site $input_line[1] for $names[$o], error: $error\n";
	    print OUTPUT "\t$alleles[0]|$alleles[1]:$end_of_line";
	} #end if
	$error=0;
    } #end for       

    $while_loop_counter++; 
} #end while                           


print OUTPUT "\n";
print STDERR "done.\n";
