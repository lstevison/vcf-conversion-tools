#! /usr/bin/perl

#program thins sites in VCF file that are too close together
#August 15th, 2012

$vcf = $ARGV[0];
$output = $ARGV[1];

unless ($#ARGV==1) {
    print STDERR "Please provide name of input and output vcf files on command line\n\n";
    die;
} #end unless

open(VCF, $vcf);
open(OUTPUT, ">$output");
$remove_both_sites = 0;
print STDERR "Reading in VCF file...";

while(<VCF>) {
    chomp;
    if ($_=~/^\#/) {
	$first_line=1;
	print OUTPUT "$_\n";
	next;
    } #end if
    
    if ($first_line==1) {
	@input_line = split(/\s+/, $_);
	$last_pos = $input_line[1];
	$last_site = $_;
	$first_line=0;
	next;
    } #end if

    @input_line = split(/\s+/, $_);
    $current_pos = $input_line[1];

    $difference = $current_pos - $last_pos + 1;

    if ($difference>15 && $remove_both_sites==0) {
	print OUTPUT "$last_site\n";
	$remove_both_sites = 0;
    } elsif ($difference>15 && $remove_both_sites==1) {
#	print STDERR "$last_pos";
	$remove_both_sites = 0;
    } elsif ($difference<=15 && $remove_both_sites==0) {
	$remove_both_sites = 1;
#	print STDERR "\n$last_pos, ";
    } elsif ($difference<=15 && $remove_both_sites==1) {
#	print STDERR "$last_pos, ";
    } #end elsif

    $last_pos = $current_pos;
    $last_site = $_;

} #end while

if ($difference>15 && $remove_both_sites==0) {
    print OUTPUT "$last_site\n";
} else {
#    print STDERR "$last_pos\n";
} #end if

print STDERR "done.\n";
