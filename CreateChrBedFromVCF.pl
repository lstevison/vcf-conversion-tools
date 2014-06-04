#! /usr/bin/perl

#Created: December 6, 2011
#Last Modified: January 7, 2013

#program reads in a vcf file and extracts coordinates
#in bed format
#of a specific number of SNPs in a window

$vcf = $ARGV[0];	
$block = $ARGV[1];

unless ($#ARGV==1) {
    print STDERR "Error: Please provide filename of input VCF file and synteny block number on command line!\n\n";
    die;
} #end unless

open(VCF, $vcf);

$SNP_counter = 0;
@chromosome = ();
@coordinates = ();
$first_time_in_loop = 1;
$first_SNP = 1;
$second_SNP = 0;

while (<VCF>) {
    chomp;
    $vcf_line = $_;
    @variant = split(/\s+/, $vcf_line);

    if ($vcf_line=~/^\#/) {
	next;			      #do not read in any header lines
    } elsif ($vcf_line=~/^chr/ && $variant[6]=~/PASS/) {
	$SNP_counter++;
    } else {
	print STDERR "Non-pass sites in VCF at $variant[0] $variant[1] bc site is $variant[6]!\n";
	next;
    } #end else
	
    if($first_time_in_loop == 1) {	 #first time in loop opens first output
	$chromosome=$variant[0];
	$first_time_in_loop = 0;
	open(OUTPUT, ">>$variant[0]-input.BED");
	open(OUTPUT2, ">$variant[0].synteny_block.$block.blocks.txt");
    } #end if
	
    if ($chromosome !~/$variant[0]/) {	 #when new chromosome, open new output
	print STDERR "finished $chromosome, now starting $variant[0]...\n";
	print OUTPUT "$last_position\n";
	print OUTPUT2 "$last_position\n";
	open(OUTPUT, ">>$variant[0]-input.BED");
	open(OUTPUT2, ">$variant[0].synteny_block.$block.blocks.txt");
	$first_SNP = 1;
	$SNP_counter=1;
    } #end if	
	
    if ($SNP_counter == 1 && $first_SNP == 1) {	 #prints first line of BED file
	print OUTPUT "$variant[0]\t$variant[1]\t" ;
	print OUTPUT2 "$variant[0]\t$variant[1]\t" ;
#	print STDERR "$variant[0]\t$variant[1]\t" ;
	$first_SNP = 0;
	$first_interval = 1;
	$after_first_interval = 0;
    } elsif ($SNP_counter == 199 && $made_past_200 ==1 && $after_first_interval==0) {				#prints end position and next start position
	print OUTPUT "$variant[1]\n$chr_next\t$start_next\t" ;
	print OUTPUT2 "$variant[1]\n$chr_next\t$start_next\t" ;
#	print STDERR "$variant[1]\n$chr_next\t$start_next\t" ;
	$made_past_200 = 0;
	$second_SNP = 1;
	$after_first_interval = 1;
    } elsif ($SNP_counter == 199 && $made_past_200 ==1 && $second_SNP==1) {				#prints end position and next start position
	print OUTPUT "$variant[1]\n$chr_next\t$start_next\t" ;
	print OUTPUT2 "$variant[1]\n$chr_next\t$start_next\t" ;
#	print STDERR "$variant[1]\n$chr_next\t$start_next\t" ;
	$made_past_200 = 0;
	$second_SNP = 0;
    } elsif ($SNP_counter == 199 && $made_past_200 ==1 && $second_SNP==0) {				#prints end position and next start position
	print OUTPUT "$variant[1]\n$chr_next\t$start_next\t" ;
	print OUTPUT2 "$variant[1]\n$chr_next\t$start_next\t" ;
#	print STDERR "$variant[1]\n$chr_next\t$start_next\t" ;
	$made_past_200 = 0;
    } elsif ($SNP_counter == 3800 && $after_first_interval == 1 && $made_it_past_200==0) {			#holds next beginning
	$chr_next = $variant[0];
	$start_next = $variant[1];
	$SNP_counter=0;
	$made_past_200 =1;
    } elsif ($SNP_counter == 3801 && $first_interval == 1) {				#holds next beginning
	$chr_next = $variant[0];
	$start_next = $variant[1];
	$SNP_counter=0;
	$made_past_200 =1;
	$first_interval = 0;
    } #end if
	
    $chromosome=$variant[0];
    $last_position = $variant[1];
#    print STDERR "Snp counter: $SNP_counter\n";    
} #end while

print OUTPUT "$last_position\n";
print OUTPUT2 "$last_position\n";
