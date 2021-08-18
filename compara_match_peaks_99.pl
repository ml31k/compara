#!/usr/bin/env perl
##### edit for general usage (202001) Huayun
# This script now connects to the public ensembl database. To connect to archived versions, the corresponding version of the API should be used.
# Edit the Registry to connect only to specific databases. Or use a text file with registry information. 
# The final output is then processed by process_compara.R
#

#### Run as perl compara_match_peaks.pl -s human -i ${path}peakfile.gff -o ${path}peakfile_out.txt -a aligntype -r species1,species2,speciesN
# ***** NOTES ***** #
#
# Notes: In order to use the ensembl modules, you need to install bioperl as well 
# For Windows Users, install Cygwin, with appropriate packages and then install bioperl and ensembl modules
# See: http://www.bioperl.org/wiki/Installing_Bioperl_on_Windows
#
# We are now always connecting to a locally installed version of ensembl on hpf
# and for the multiple alignments we are using either "EPO" or "PECAN"
# 
# If you wish to use a specific set of alignments from a specific ensEMBL version
# you can use the "getComparaDBIDs.pl" script.
# 
# The final output is to be used by process_compara.R
#
# ***** END ***** #


#use Modern::Perl 2011;
use autodie;

# Use the current perl api version 
# Target libraries included in lib folder. Set path accordingly

my $comparaPath = "/hpf/largeprojects/mdwilson/scripts/Compara/scripts";

#use lib "$comparaPath/lib/compara99/ensembl/modules/";
#use lib "$comparaPath/lib/compara99/ensembl-compara/modules/";
#use lib "$comparaPath/lib/";

use lib "/hpf/largeprojects/mdwilson/scripts/Compara/scripts/lib/";;
use lib "/hpf/largeprojects/mdwilson/scripts/Compara/scripts/lib/compara99/ensembl/modules/";;
use lib "/hpf/largeprojects/mdwilson/scripts/Compara/scripts/lib/compara99/ensembl-compara/modules/";;

use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::SimpleAlign;
use Bio::AlignIO;
use Bio::LocatableSeq;
use Getopt::Long;
use Data::Dumper;

use List::MoreUtils;
use Data::Dump qw(pp); 
pp(\%INC);


sub trim($) {
	my $string =shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}


# If something goes wrong then the original command call will be returned to STDERR
warn join(" ","compara_match_peaks.pl",@ARGV),"\n";
my $helpmessage= "perl compara_match_peaks.pl -s species -i infile -o outfile -q seqfile -r species1,species2,speciesN -a EPO|PECAN \nList of species: human, macaque, chimp, orangutan, gorilla, marmoset, mouse, rat, dog, chicken \nensembl compara release 99\n";
 
&GetOptions(
    'outfile|o:s'              => \$outfile,
    'seqfile|q:s'              => \$seqfile, # where to store MSA sequences
    'species|s:s'              => \$species, # peak file species
    'infile|i:s'               => \$infile,
    'alignment_type|a:s'	   =>\$alignment_type,
    'species_out|r:s'		   => \$species_out, #common names of species to be used in output
    'help|?'                   => \$help,
    );

## Michael 2014.2.04:
# Wrapper will input appropriate db_names directly. Manual run of this script should likewise use db_names.
#Replaced non-critical 'throw' statements with warn - next() for better stability

my @dbnames_used = split("\\,",$species_out);

die "USAGE: $helpmessage" if $help;
die "USAGE: $helpmessage"  if ((not defined $outfile) ||(not defined $infile)||(not defined $species)||(not defined $species_out)||(not defined $alignment_type)||(not defined $seqfile)) ;

open(IN, "$infile") or die("Cannot open $infile");
open(OUT, ">$outfile") or die("Cannot write $outfile");
open(SEQ, ">$seqfile") or die("Cannot write $seqfile"); 

## Print header in outfile

print OUT join("\t","#Original_peak" , @dbnames_used)."\n";

#-- This is the minimum number of BP in an alignment
my $min_bp_in_almnt = 10;
my $coord_system = "chromosome";


#########################################################
# Registry
# Bio::EnsEMBL::Registry->load_registry_from_db(
#     -host => 'ensembldb.ensembl.org', 
#     -user => 'anonymous'
#     );

#my $registry = "Bio::EnsEMBL::Registry";
Bio::EnsEMBL::Registry->load_registry_from_db(
    -host => 'ensembldb.ccm.sickkids.ca', 
    -user => 'readonly',
    -pass => 'readonly',
    -verbose => "1"
    );
# Alternative: registry with a txt file with registry information
# my $registry="ebi_ensg67_compara_registry.txt"; 
# Bio::EnsEMBL::Registry->load_all($registry);

#########################################################
# Specify the kind of multi-species alignment we want to use
# Currently using "13 eutherian mammals (EPO)" and PECAN

$alignment_type = "EPO" unless ($alignment_type);
if($alignment_type eq "EPO"){
    $set_of_species = "mammals"
} elsif ($alignment_type eq "PECAN") {
    $set_of_species = "amniotes"
};

# Get MethodLinkSpeciesSet adaptor. This adaptor links various analyses (method_link_type) with a set of species (species_set)
my $method_link_species_set_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi', 'compara', 'MethodLinkSpeciesSet');

throw("Cannot connect to Compara") if (!$method_link_species_set_adaptor);

# Get the method_link_species_set for multiple alignments of mammals
#my $set_of_species = $method_link_species_set_adaptor->species_set();
#print $set_of_species;

my $method_link_species_set = $method_link_species_set_adaptor->fetch_by_method_link_type_species_set_name($alignment_type, $set_of_species);

throw("The database does not contain any $alignment_type data for $set_of_species!")
    if (!$method_link_species_set);
print "I got the mlss for ", $method_link_species_set->name(),"! and it covers the following species: \n";

# For the method:species_set combination, get the genome information of all the species 
my $mlss_gdbs = $method_link_species_set->species_set->genome_dbs();
# Output the names of all species in the species_set
foreach $mls (@{$mlss_gdbs}) {
  print $mls->name(), "\n";
}

#########################################################
# Getting the adaptors for the core database of the species and genomic align blocks

my $slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'Slice');
throw("Registry configuration file has no data for connecting to <$species>") if (!$slice_adaptor);
print "All slice adaptors loaded!\n";

my $genomic_align_block_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi', 'compara', 'GenomicAlignBlock');
throw("Cannot get GenomicAlignBlock adaptors") if (!$genomic_align_block_adaptor);
print "All genomic align blocks loaded!\n";

my $start_time = localtime();
print "Starting at $start_time"."\n";
#########################################################
#loop through the peak file, get slice and get alignment block for each slice, then restrict between alignment positions #
while (my $peak = <IN>) {
    # remove erroneous line break 
    chomp($peak);
    # remove white space using self-defined trim subroutine
    trim($peak);
    # Split the "peak" into an array by tabs
    my @peak = split(/\t/,$peak);
    
    ## Recover peak information
    my $peak_chr = $peak[0];
    my $peak_start = $peak[1];
    my $peak_end = $peak[2];
    
	# Generate query slice
    my $query_slice = $slice_adaptor->fetch_by_region('toplevel', $peak_chr, $peak_start, $peak_end);
    my $query_peak = $species.".".$peak_chr.".".$peak_start.".".$peak_end;

    if (!$query_slice)
    {
	warn("No Slice can be created with coordinates $peak_chr:$peak_start-". "$peak_end");
	next();
    };

	# Initiate hash for MSA
    my %orthologous_regions=(); # To store MSA regions
    my %orthologous_sequences=(); # To store MSA sequences
    
    if ($query_slice) {
	# query the GenomicAlignBlock adaptor using the slice we retrieved for the peak. 
	
	# IMPORTANT NOTE 120816
	# This call can generate an array of genomic align blocks of size greater than one as well!
	# As a result each genomic_align_blocks is an array of references to genomic_align_arrays
	# if the size of $genomic_align_blocks > 1 then there are multiple different aligned segments
	# and these will end up generating multiple compara peak data points! BEWARE!!!

	my $genomic_align_blocks = $genomic_align_block_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($method_link_species_set, $query_slice);
	if (!$genomic_align_blocks)
	{
	    warn("Cannot get GenomicAlignBlock");
	    next();
	};
 
	##### IMPORTANT ADDITION #####
	# 120813: Calvin Mok
	# Sometimes the Compara database will return NOTHING even though a specific slice is valid
	# In previous iterations of the program, it would result in a discrepancy between
	# the final number of "Compara peaks" vs the original set of .gff peak data
	# In these cases, I have now included the original peak just in case
	# This may skew the data so it can be removed as well  
	my $hits = 0;
	
	#Print the original peak.
	print OUT $query_peak."\t";
	
	if (scalar @{$genomic_align_blocks} < 1) {
	    $hits = 1;
	}
	
	
	foreach my $this_genomic_align_block (@$genomic_align_blocks) {
	    # Correcting the position of genomic align block: restrict to initial peak region
	    $this_genomic_align_block_aux= $this_genomic_align_block->restrict_between_reference_positions($peak_start,$peak_end);
	    
	    if (!$this_genomic_align_block_aux){
		warn ("Could not restrict between reference position for $peak_chr: $peak_start to $peak_end \n") ;
		next();
	    }
	    else{
		$this_genomic_align_block=$this_genomic_align_block_aux;
	    }

		# Get all the genomic alignments within the genome align block
	    if ($this_genomic_align_block) {
		my $genomic_align_array = $this_genomic_align_block->genomic_align_array();
		#print Dumper $genomic_align_array; 
		
		foreach my $genomic_align (@$genomic_align_array) {
		    if ((List::MoreUtils::first_index { $_ eq $genomic_align->genome_db()->name() } @dbnames_used) != -1) {
			# loops through for each species that is pulled down for that fragment and only use the alignments for specified species
			
			$hits += 1;
			## This was added to clean the genome fragment name, since for rat sometimes includes the version
			## ID AABR06109730 amedina20130108
			## Is confunsing to have things mapping to the external contigs, if this match to a weird chromosome is going to be dropped (next) amedina 20130110
		     
			my $genome_id=$genomic_align->genome_db()->name();
			my $genome_fragment=$genomic_align->dnafrag->name();
			if($genome_fragment=~/\./){
			    $genome_fragment =~s/\..*//;
##			    next;
			}
		
			# Get the coordinate of the aligned region and save to hash	
			my $align_coord = $genome_id.".". $genome_fragment.".".$genomic_align->dnafrag_start().".".$genomic_align->dnafrag_end(); 
			push(@{$orthologous_regions{$genome_id}},$align_coord );
			
			# ML:2015.06.16
			# Added a feature to retrieve the actual MSA sequences for each peak/fragment that is present. May be useful for downstream applications. If it's a deletion (returns a region of -1 length, skip)
			unless ($orthologous_sequences{$query_peak}->{$genome_id} || $genomic_align->dnafrag_start > $genomic_align->dnafrag_end)
			{
			    $orthologous_sequences{$query_peak}->{$genome_id} = "$align_coord:".$genomic_align->aligned_sequence();
			}
			 
		    }
		}
	    }
	}
	# Output the results. Add tabs for importing into R with stable formatting
	my $ortho_regions="";
	my $ortho_seqs="";

	foreach my $genomeID ( @dbnames_used){
	    if ( $orthologous_regions{$genomeID}){
		$ortho_regions.=join(":", @{$orthologous_regions{$genomeID}})."\t";
	    }
	    else {
		$ortho_regions.="NA"."\t";
	    }
	}
	$ortho_regions=~s/\t$//;
	print OUT $ortho_regions."\n";


	foreach my $seq_peak (keys %orthologous_sequences){
	    $ortho_seqs.="\n>$seq_peak\n";
	    foreach my $genomeID ( @dbnames_used){
		if ($orthologous_sequences{$seq_peak}->{$genomeID}) {
		    $ortho_seqs.=$orthologous_sequences{$seq_peak}->{$genomeID}."\n";
		} else {
		    $ortho_seqs.=$genomeID.".NA.NA.NA:\n";
		}
	    }
	}
	print SEQ $ortho_seqs;
    }
}

close IN;
close OUT;

my $end_time = localtime();
print "Ending at $end_time"."\n";

print "Results stored in ".$outfile."\n";
print "MSA seq stored in ".$seqfile."\n";
exit;

