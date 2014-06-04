#! /usr/bin/perl

#Program reads in MS formatted input (*.hud) and converts to LDhat input format

$input=$ARGV[0];
$output=$ARGV[1];
$num_ind=$ARGV[2];
$num_haps=$num_ind*2;

unless ($#ARGV==2) {
    print STDERR "Please provide input and output filenames and number of individuals on command line.\n\n";
    die;
} #end unless 

open(INPUT, $input);
open(SITES, ">$output.ldhat.sites");
open(LOCS, ">$output.ldhat.locs");

$num_sites=(<INPUT>);
chomp $num_sites;
$positions=(<INPUT>);
chomp $positions;
@positions=split(/\s+/,$positions);

print LOCS "$num_sites\t$positions[$#positions]\tL\n";

for ($p=0; $p<=$#positions;$p++) {
    $kb_sites = $positions[$p]/1000;
    printf LOCS ("%0.3f\n", $kb_sites);
#    print LOCS "$positions[$p]\n";
} #end for

#print STDERR "Num sites: $num_sites; Length positions array: $#positions.\n";
print SITES "$num_haps\t$num_sites\t1\n";

$loop_counter=1;

while (<INPUT>) {
    chomp;
    if ($loop_counter >$num_haps) {
	last;
    } #end if

    print SITES ">haplotype $loop_counter\n$_\n";
    $loop_counter++;

} #end while
