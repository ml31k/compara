#!/hpf/tools/centos6/perl/5.18.2/bin/perl
##### edit for hpf (121022) 
# fixed output format (for R read.table)
# alignment is now an option (EPO/PECAN)
# compara species -r option (indexing species names)


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

# Current perl api version 70
use lib "/hpf/projects/mdwilson/lib/comparaAPI/compara70/ensembl/modules";
use lib "/hpf/projects/mdwilson/lib/comparaAPI/compara70/ensembl-compara/modules";
use lib "/hpf/projects/mdwilson/lib/bioperl-1.2.3";

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

sub trim($) {
	my $string =shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}


# If something goes wrong then the original command call will be returned to STDERR
warn join(" ","compara_match_peaks.pl",@ARGV),"\n";
my $helpmessage= "perl compara_match_peaks.pl -s species -i infile -o outfile -q seqfile -r species1,species2,speciesN -a EPO|PECAN \nList of species: human, macaque, chimp, orangutan, gorilla, marmoset, mouse, rat, dog, chicken \nensembl compara release 70\n";
 
&GetOptions(
    'outfile|o:s'              => \$outfile,
    'seqfile|q:s'              => \$seqfile, # where to store MSA sequences
    'species|s:s'              => \$species, # peak file species
    'infile|i:s'              => \$infile,
    'alignment_type|a:s'			=>\$alignment_type,
    'species_out|r:s'			=> \$species_out, #common names of species to be used in output
    'help|?' => \$help,
    );
################
#Sepecies parse common names
#dog name in ensembl is canfam
#$species=~ s/dog/canfam/;
#$species_out=~ s/dog/canfam/;


### species info is generated for the compara species selected in -r/-s options
# Species information (latin/common/db naming schemes)
# Subsets of these files are generated to include only those wanted in the analysis (-r option)

## read table from /home/mattc/scripts/OrgNomenclature.txt to replace lines below
#my @commonnames = ("human","macaque","chimp","orangutan","gorilla","marmoset","mouse","rat","canfam","chicken");		
#my @latinnames = ('Homo sapiens','Macaca mulatta','Pan troglodytes','Pongo pygmaeus','Gorilla gorilla','Callithrix jacchus','Mus musculus','Rattus norvegicus','Canis familiaris','Gallus gallus');
#my @dbnames = ('homo_sapiens','macaca_mulatta','pan_troglodytes','pongo_pygmaeus','gorilla_gorilla','callithrix_jacchus','mus_musculus','rattus_norvegicus','canis_familiaris','gallus_gallus');

## return index of peak file species (from -s option)
#my $peakspec_index = List::MoreUtils::first_index { $_ eq $species } @commonnames;

## return index of species selected in compara (-r option)
#my @species_used = split("\\,",$species_out);
#my @used_index = (1..scalar(@species_used));#initialize used_index
#foreach $i (0..(scalar(@species_used)-1)){
#	@used_index[$i] =  List::MoreUtils::first_index { $_ eq $species_used[$i] } @commonnames;
#};

##create subset of naming tables for species in -r option
#my @commonnames_used=@commonnames[@used_index];
#my @latinnames_used = @latinnames[@used_index];
#my @dbnames_used = @dbnames[@used_index];
#
#print "Species peak file: " . $species."\n";
#print "Species to be reported: ". join ("\t",@commonnames_used)."\n\t".join("\t", @latinnames_used)."\n\t".join("\t",@dbnames_used )."\n";

## Michael 2014.2.04:
#Removed common name parsing to run with a wrapper. Wrapper will input appropriate db_names directly. Manual run of this script should likewise use db_names.
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

# necessary if code under does not work where ebi_ensg67_compara_registry is a txt file with registry information
# my $registry="ebi_ensg67_compara_registry.txt"; 
################################################################
#####   check this file exists in the installation of the local ensembl in the web is called registry configuration file 
# Bio::EnsEMBL::Registry->load_all($registry);

Bio::EnsEMBL::Registry->load_registry_from_db(
    -host => 'mysqld1', 
    -user => 'ensemdba',
    -port => 3311,
	-pass => 'ensemDB@',
    -db_version => 70,
    -verbose => 1
    );

# Here is where we specify the KIND of multi-species alignment we want to make
# In order to add alignments you can find the align_dbID (alignment ID number) you wish to use with getComparaDBIDs.pl
# and the align_species_set from http://www.ensembl.org/info/docs/compara/analyses.html
# currently using "13 eutherian mammals (EPO)" and PECAN

my %align_species_set = ("EPO" => 'mammals',"PECAN" => 'amniotes');		
my %align_dbID = ("EPO" => 619,"PECAN" => 620);
$alignment_type="EPO" unless ($alignment_type);
my $dbID = $align_dbID{$alignment_type};
my $species_set_name= $align_species_set{$alignment_type};


my $method_link_species_set_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi', 'compara', 'MethodLinkSpeciesSet'); #my database, multi using compare, using method to link species

#print "$method_link_species_set_adaptor\n";
#print "$alignment_type is $dbID\n";

my $method_link_species_set = $method_link_species_set_adaptor->fetch_by_dbID($dbID);

throw("The database does not contain any $alignment_type data for $set_of_species!")
    if (!$method_link_species_set);
print "I got the mlss for ", $method_link_species_set->name(),"! and it covers the following species: \n";

my $mlss_gdbs = $method_link_species_set->species_set_obj->genome_dbs();

foreach $mls (@{$mlss_gdbs}) {
  print $mls->name(), "\n";
}


#### Getting the adaptors
# In this case, the $species value would normally be 'human' from command line initiation
# I'm guessing later that you'll only do direct human/other comparisons through the program
# And not a cross-species alignment
my $slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'Slice');
throw("Registry configuration file has no data for connecting to <$species>") if (!$slice_adaptor);
print "All slice adaptors loaded!\n";

my $genomic_align_block_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi', 'compara', 'GenomicAlignBlock');
#print $genomic_align_block_adaptor; die "BOOM";
throw("Cannot get GenomicAlignBlock adaptors") if (!$genomic_align_block_adaptor);
print "All genomic align blocks loaded!\n";

#loop through the peak file, get slice and get alignment block for each slice, then restrict between alignment positions #
while (my $peak=<IN>) {
    # remove erroneous line break 
    chomp($peak);
    # remove white space using self-defined trim subroutine
    trim($peak);
    # Split the "peak" into an array by tabs
    my @peak=split(/\t/,$peak);
    
    ## Recover peak information
    my $peak_chr=$peak[0];
    my $peak_start=  $peak[1];
    my $peak_end= $peak[2];
    
    my $query_slice = $slice_adaptor->fetch_by_region('toplevel', $peak_chr, $peak_start, $peak_end);
    my $query_peak = $species.".".$peak_chr.".".$peak_start.".".$peak_end;
#    print "grabbed query slice $peak_chr $peak_start $peak_end\n" ;
    if (!$query_slice)
    {
	warn("No Slice can be created with coordinates $peak_chr:$peak_start-". "$peak_end");
	next();
    };
    my %orthologous_regions=(); #Added to store MSA sequences
    my %orthologous_sequences=();
    
    if ($query_slice) {
	# print "definitely a query slice $peak_chr $peak_start $peak_end\n" ;
	# query the adaptor using the slice we retrieved for the peak. Not sure if this exists in GenomicAlignBlock
	# this code is nearly identical to the compara tutorial code so it should be correct
	
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

#	print "size of alignment block: ", scalar @{$genomic_align_blocks},"\n";    
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
	    
	    #print Dumper $this_genomic_align_block;
	    #die "+".$this_genomic_align_block."+";
	    my $genomic_align_array = $this_genomic_align_block->genomic_align_array();
#	    print "size of alignment array: ", scalar @{$genomic_align_array},"\n";
	    #print Dumper $genomic_align_array;
	    #die;
	    
	    ################
	    ## Code to debbug alignment. It loops thorugh each line of the aligment, printing the organism name.
	    #foreach my $genomic_align (@$genomic_align_array) {
	    # sanity check to see if the genomic alignments match with our defined species
	    # print ("Sanity check: ", $genomic_align->genome_db()->name(), " vs. ", $dbnames[$peakspec_index], "\n");
	    #print Dumper  $genomic_align;
	    #if (($genomic_align->genome_db()->name()) eq ($dbnames[$peakspec_index])) {
	    # print ("accepted\n");
	    # print Dumper $species_xtra{$genomic_align->genome_db()->name()};
	    # print Dumper $species_ok{$species};
	    # print "FOUND\n";
	    # print Dumper $genomic_align->genome_db()->name().".".$genomic_align->dnafrag->name().".".$genomic_align->dnafrag_start(),"\n";
	    
	    # exits looping through the genomic align blocks if it finds a human vs. human alignment
	    #last;
	    # }
	    #die;
	    #}
	    #die;
	    ################


	    #Correcting the position: restrict to initial peak region
	    $this_genomic_align_block_aux= $this_genomic_align_block->restrict_between_reference_positions($peak_start,$peak_end);
	    
	    if (!$this_genomic_align_block_aux){
		warn ("Could not restrict between reference position for $peak_chr: $peak_start to $peak_end \n") ;
		next();
	    }
	    else{
		$this_genomic_align_block=$this_genomic_align_block_aux;
	    }
	    if ($this_genomic_align_block) {
		my $genomic_align_array = $this_genomic_align_block->genomic_align_array();
		#print Dumper $genomic_align_array; 
		
		foreach my $genomic_align (@$genomic_align_array) {
		    #print Dumper $genomic_align->genome_db()->name();
		    # If the alignment isn't in our species list, we dump it
		    #print ("Final sanity check for species: ", $genomic_align->genome_db()->name(), "\n");
		    if ((List::MoreUtils::first_index { $_ eq $genomic_align->genome_db()->name() } @dbnames_used) != -1) {
			# output looks like a genome dbname.alignment fragment name.start.end
			# loops through for each species that is pulled down for that fragment
			
			$hits += 1;
			## This was added to clean the genome fragment name, since for rat somtimes includes the version
			## ID AABR06109730 amedina20130108
			## Is confunsing to have things mapping to the external contigs, if this match to a weird chromosome
			## is going to be dropped (next) amedina 20130110
		     
			my $genome_id=$genomic_align->genome_db()->name();
			my $genome_fragment=$genomic_align->dnafrag->name();
			if($genome_fragment=~/\./){
			    #my $genome_fragment2=$genome_fragment;
			    #$genome_fragment2=~s/\./_/; # contains the name of the organims
			    #$genome_fragment=$genome_fragment2;
			    #print $genome_id.".". $genome_fragment.".".$genomic_align->dnafrag_start().".".$genomic_align->dnafrag_end(); 
			    #die "BOOM";
			    next;
			}
			##
			
			my $align_coord = $genome_id.".". $genome_fragment.".".$genomic_align->dnafrag_start().".".$genomic_align->dnafrag_end(); 

			# ML:2015.06.16
			# Added a feature to retrieve the actual MSA sequences for each peak/fragment that is present. May be useful for downstream applications.
			unless ($orthologous_sequences{$query_peak}->{$genome_id})
			{
			    $orthologous_sequences{$query_peak}->{$genome_id} = "$align_coord:".$genomic_align->aligned_sequence();
			}
			push(@{$orthologous_regions{$genome_id}},$align_coord );

			#print $align_coord."--\t---";
			
			#print OUT $align_coord."\t";
		           
			
		     
		    }
		}
	    }
	}
	# Adding tabs for importing into R with stable formatting
	# Errors were caused in earlier versions since R did not know the correct number of columns from this file
	#print OUT "\t" x (scalar(@species_used)-$hits-1);
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

print "Results stored in ".$outfile."\n";
print "MSA seq stored in ".$seqfile."\n";
exit;


## print OUT $dbnames{$species}.".".$peak_chr.".".$peak_start.".".$peak_end;

##  	    	if ($species_xtra{$genomic_align->genome_db()->name()} eq $species_ok{$species}) {

#  	       if ($species_xtra{$genomic_align->genome_db()->name()}) {

#	print OUT "\t" x (scalar(keys %species)-$hits-1);
