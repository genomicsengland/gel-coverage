#this script generates bed file with coordinates of DNA fragments that do not contain Ns
#only autosomes are analysed
#similar file for GRCh37 is at /accelrys/apps/gel/genomeref/data/human/human_g1k_v37_NonN_Regions.CHR.bed
#Razvan complimentary file is at /genomes/resources/genomeref/data/human/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.N_regions.bed

$f = 0;
$chrom = "";
#$ref = "/genomes/resources/genomeref/Illumina/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa";
#$bed = "/genomes/scratch/asosinsky/resources/GRCh38Decoy_NonN_Regions.CHR.bed";

$ref = "/genomes/resources/genomeref/data/human/GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa";
$bed = "/genomes/scratch/asosinsky/resources/Homo_sapiens.GRCh37.75.dna.primary_assembly.NonN_Regions.CHR.bed";

open(REF, "$ref");
open(BED, ">$bed");
while($line = <REF>){
	chomp($line);
	#if ($line =~ /^>(chr\d+)\s+(.*)/){
	if ($line =~ /^>(\d+)\s+(.*)/){
		deploy($chrom, $start, $coord) if $f && $chrom;
		$chrom = $1;
		$coord = 0;
	}elsif($line =~ /^>(\S+)\s+(.*)/){
		deploy($chrom, $start, $coord) if $f && $chrom;	
		$chrom = "";
		$coord = 0;
	}elsif($line =~ /^N+$/i && $chrom){
		deploy($chrom, $start, $coord) if $f && $chrom;
		$coord += length($line);
	}elsif($line =~ /^([TACGMKYRSWVBDHN]+)$/i && $chrom){
		$seq = $1;
		@seq = split(//,$seq);
		foreach $nucl (@seq){
			if ($nucl !~ /N/i && $f == 0){
				$start = $coord;
				$f = 1;
			}elsif ($nucl  =~ /N/i){
				deploy($chrom, $start, $coord) if $f && $chrom;
			}
			$coord ++;
		}

	}elsif($line =~ /^([TACGMKYRSWVBDHNtacgmkyrswvbdhn]+)$/){
		
	}else{
		print "Unrecognised format: $line\n";
	}
}
deploy($chrom, $start, $coord) if $f && $chrom;

sub deploy {
	($chrom, $start, $coord) = @_;
	print BED "$chrom\t$start\t$coord\n" ; 
	$f = 0; $start = 0;
}