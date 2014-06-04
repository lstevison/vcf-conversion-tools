vcf-conversion-tools
====================

Tools to (1) convert to and from vcf format, (2) convert between LDhat, fastPHASE, and PHASE formats, and (3) convert VCF to BED format and get overlapping bins of 4000 or 400 SNPs.

(1) Files in this folder will help convert between vcf and other formats:
'fastPHASE2VCF.pl'  #converts between fastPHASE output to VCF
'thinVCF.pl'   #thins VCF files (slightly different algorithm from vcftools, removes all but one site if sites are close to each other)
'vcf_merge.pl'  #merges multiple VCF files into single file 
'vcf2fastPHASE.pl'  #converts a VCF file to fastPHASE input for autosomes and females
'vcf2fastPHASE_4males.pl'   #converts VCF file to fastPHASE input for males on the X
'vcf2MS.pl'   #converts VCF to MS format


(2) A few other programs are included as they an be used with the others to convert to and from LDhat formats: 
'fastPHASE2LDhat.pl'  #converts fastPHASE output to LDhat input
'MS2LDhat.pl'   #converts MS input to LDhat input
'MS2PHASE.pl'   #converts MS input to PHASE input

(3) The following programs are useful for creating bed files of a subset of overlapping SNPs from a VCF file:
'CreateChrBedFromVCF.pl'  #creates BED file with bins of 4k SNPs, 100 overlapping, useful for running LDhat
'CreateChrBedFromVCF2.pl'   #creates BED file with bins of 400 SNPs, 100 overlapping, useful for running PHASE