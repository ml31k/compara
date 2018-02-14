#!/usr/bin/perl

## Michael - 2014.04.04 V.1.01 Changes:
## Changed outspecies option from -o to -out_spec; added -o outputdirectory

## This script is designed to function as a wrapper for the compara package.
## It will build the directory structure and appropriately name files for
## subsequent analysis.
## It will generate text files with commands that can subsequently be submitted
## to the cluster

use strict;
use Getopt::Long;
use Cwd;
use Cwd 'abs_path';
use File::Path 'make_path';

my $cwd = getcwd();
my $outdir = "$cwd";
my $inputFile;
my $help;
my $path = "/hpf/largeprojects/mdwilson/scripts/Compara/scripts";
my $out_species;
my %species_names;
my $is_gff;
my $alignment = "EPO";
my $mmr = "0:1";
my $do_paralogs = "FALSE";
my $auto_queue;
my $command_prefix = " ";
my $out_species_string;
my $priority_species;
my @output_species;

GetOptions
    ("h|help" => \$help,
     "i=s" => \$inputFile,
     "o=s" => \$outdir,
     "out_spec=s" => \$out_species,
#     "is_gff" => \$is_gff,
     "a=s" => \$alignment,
     "mmr=s" =>\$mmr,
     "p=s" => \$do_paralogs,
     "q=s" => \$auto_queue,
     "p_spec=s" => \$priority_species,
     "cmd=s" => \$command_prefix
     )
    or die ("invalid commandline args\n");

my $helpmsg = join ("\n", 
		    "##############################################",
		    "Compara wrapper version 1.01 by Mimi31k",
		    "##############################################",
		    "Usage: -i inputfile [OPTIONS]",
		    "",
		    "Mandatory Arguments:",
		    "\t-i input\t\tA tab delimited file with columns Species ID(ex. TF name) path/to/peak",
		    "\t-o outdir\t\tABSOLUTE PATH to Output directory. Defaults to CURRENT_WORKING_DIRECTORY",
		    "\t\t\t\tSpecies names must correspond to species_db.txt found in the compara directory",
		    "\t\t\t\tPeak files may be in any format with first 3 columns chr# Start End (no header)",
		    "\t\t\t\tNon .gff files will be converted and stored locally. gff files will be linked.",
		    "",
		    "Optional Arguments:",
		    "\t-h | -help\t\tOpen this help menu and exit",
#		    "\t-is_gff\t\t\tPeak files are in .gff format and do not require conversion. False by default",
		    "\t-a alignment\t\tAlignment mode to run Compara [EPO|PECAN]. EPO by default",
		    "\t-mmr max:min\t\tMax:Min ratios for peak overlaps. Default 0:1",
		    "\t-p do_paralogs\t\t[TRUE] Determines if R component runs paralog analysis. FALSE by default",
		    "\t-q auto_queue\t\t[FALSE] Determines if compara perl commands are sent to qsub. TRUE by default",
		    "\t-cmd command prefix\tCommand prefix for all compara commands (essential source commands, etc.)",
		    "\t-out_spec out_species\t\tComma delimited list of output species. Equal to input species by default",
		    "\t-p_spec priority_species\t\tName of species that should always comes first in Compara output",
		    "",
		    "Output_species_order: Custom order spec. by -out_spec >",
		    "Priority species first + alphabetized > alphabetized input species order (Default)\n"
    );

if ($help || (not defined $inputFile) || ($alignment && $alignment != "EPO" && $alignment ne "PECAN") || ($do_paralogs ne "FALSE" && $do_paralogs ne "TRUE") || ($auto_queue && $auto_queue ne "FALSE"))
{
    print $helpmsg;
    exit;
};
## Auto queueing of commands is set on by default, unless specified as "FALSE"
$auto_queue = $command_prefix unless $auto_queue == "FALSE";

$outdir = abs_path($outdir);

unless ( -e $outdir )
{
    make_path($outdir);
};

## Loads the species database and common name conversion table. This file is found in
## the compara directory and can be edited to add new species/shorthand. Be sure to use
## the correct database name for compara and seperate all shorthands with a ;
open SPECIES_DB, $path."/species_db.txt" or die ("couldn't load species names!\n");
while (<SPECIES_DB>)
{
    chomp();
    my @tmp= split(" ", $_);
    $species_names{$tmp[0]} = $tmp[1];
};
close (SPECIES_DB);


## Loads the output species list and transforms them to appropriate DB names
## This is run only if an out species file is designated.
## Default is to use the same set of species as in the input file.
if ($out_species)
{
    @output_species = split(/,/, $out_species); 
    &renameSpecies (@output_species);
    $out_species_string = join(',', @output_species);
};

## Parses the input file, transforms input species to appropriate DB names
## Output species, if not already defined, are defined at this step
my %hash_by_IDs;

open INPUT, $inputFile or die ("couldn't open inputfile $inputFile");
while (<INPUT>)
{
    chomp();
    my($a_species, $a_ID, $a_peakFile) = split(" ", $_);
    &renameSpecies($a_species);    
    $out_species_string .= "$a_species," unless $out_species || $out_species_string =~ $a_species;
    if ($hash_by_IDs{$a_ID})
    {
	push (@{$hash_by_IDs{$a_ID}}, [$a_species, abs_path($a_peakFile)]);
    }
    else 
    {
	$hash_by_IDs{$a_ID} = [];
	push (@{$hash_by_IDs{$a_ID}}, [$a_species, abs_path($a_peakFile)]);
    };
};
close (INPUT);

## Determine the order of species output in Compara file:
## Priority: Custom order (as specified. by out_spec)> Priority species first > input species order (Default)

if ($out_species)
{
    print STDERR "User specified out species order: $out_species_string\n";
}
else
{
    ## Presort the list alphabetically.
    $out_species_string = join(",", sort(split(/,/, $out_species_string)));
    
    if ($priority_species)
    {
	print STDERR"User specified priority species: $priority_species\n";
	&renameSpecies($priority_species);
	print STDERR"Renamed as: $priority_species\n";
	
	## Pops out priority species and puts it at the front of list.
	$priority_species .= ",";
	$out_species_string =~ s/$priority_species//g;
	$out_species_string = $priority_species.$out_species_string;
	print STDERR "Prioritized out species + alphabetized order: $out_species_string\n";
    }
    else
    {
	print STDERR "No out species order specified. Using alphabetized input species order: $out_species_string\n";
    };
}
@output_species = split(/,/, $out_species_string);

## Setup subdirectory structure and generates peaks if necessary.
## Also prints the appropriate commands to run the perl and R portions of compara
for my $ID (keys %hash_by_IDs)
{
    if (-d "$outdir/$ID")
    {
	print STDERR "directory $ID exists!\n";
    }
    else
    {
	system ("mkdir $outdir/$ID") == 0 or die ("failed to make subdirectory for $ID\n");
	system ("mkdir $outdir/$ID/peaks") == 0 or die ("failed to make subdirectory for $ID\n");
    };
    open CMD_OUT, ">$outdir/$ID/run_compara_cmd.txt" or die("failed to write CMD_OUT commands\n");
    open R_INPUT, ">$outdir/$ID/rscript_input.txt" or die("failed to write R_INPUT commands\n");
    open R_INSTRUCT, ">$outdir/$ID/rscript_instructions.txt" or die("failed R_INSTRUCT to write commands\n");
    open R_CMD, ">$outdir/$ID/rscript_run.R" or die("failed to write R_CMD commands\n");
    my $ID_list = \@{$hash_by_IDs{$ID}};    
    my %r_input;
    
    for my $spec_and_file (@{$ID_list}) 
    {
	my $peakFileIn = $spec_and_file -> [1]; #The spec_and_file variable holds in [0] species name; [1] path/to/peakfile
	my $fileOut =  $ID."_".$spec_and_file->[0];      
	
	## Link peaks if gff, else generate gff files. named: Species_db_name_ID.gff
	if ($is_gff)
	{
	    system ("ln -s $peakFileIn $outdir/$ID/peaks/$fileOut.gff") == 0 or die ("failed to link $peakFileIn\n"); 
	}
	else
	{
	    &makeGFF ($peakFileIn, "$outdir/$ID/peaks/$fileOut.gff");
	};
	
	## Generate command to run compara:
	print CMD_OUT ("$command_prefix ",
		       "cd $outdir ; $path/compara_match_peaks.pl ",
		       "-i $outdir/$ID/peaks/$fileOut.gff ",
		       "-o $outdir/$fileOut.compara.regions ",
		       "-q $outdir/$fileOut.MSA.sequences ",
		       "-s $spec_and_file->[0] ",
		       "-r $out_species_string ",
		       "-a $alignment\n");
	## Generate the input file list for R component
	$r_input{$spec_and_file->[0]} = "$spec_and_file->[0]\t$outdir/$fileOut.compara.regions\t$outdir/$ID/peaks/$fileOut.gff\n";
    };
    ## Generate the instruction file for R script
    print R_INSTRUCT "$outdir/$ID/rscript_input.txt\t$ID\t$mmr\n";
    ## Generate the actual R script to be executed
    print R_CMD ("source (\"$path/peakConverter.R\")\n",
		 "source (\"$path/process_compara.R\")\n",
		 "source (\"$path/batchCompara.R\")\n",
		 "do.paralogs <<- $do_paralogs\n",
		 "batchComparaData(\"$outdir/$ID/rscript_instructions.txt\")\n");
    ## Prints R_INPUT in order specified by out_species_string;
    for my $out_sp (@output_species)
    {
	print R_INPUT $r_input{$out_sp};
    };
    close (R_INPUT);
    close (R_INSTRUCT);
    close (R_CMD);
    close (CMD_OUT);
    ## Submits perl commands to the cluster unless auto_queue was specified as "FALSE"
    &queue_commands("$outdir/$ID/run_compara_cmd.txt", $ID) unless $auto_queue;    
};

print "Wrapper function completed!\n";
exit;

### Link peak files or create .gffs if necessary:
#if ($is_gff)
#{
#    for (my $i; $i <= $#ID_list; $i ++)
#    {
#	system ("ln -s $peakFiles[$i] $ID[$i]/peaks/$in_species[$i]_$ID[$i].gff") == 0 or die ("failed to link $peakFiles[$i]\n");
#    };
#}
#else 
#{
#    &makeGFF ($peakFiles[$i], "$ID[$i]/peaks/$in_species[$i]_$ID[$i].gff")
#};

sub queue_commands
## Submits commands for perl compara scripts directly to the cluster. Modify this as per your submission methods
## Accepts file containing list of commands as 1st arg, ID as second.
{
    system ("queue_commands.py -i $_[0] -prefix $_[1]compara -qsub_memory 13g -q long") == 0 or die ("Error in submitting jobs to cluster!\n");
    print "System command: queue_commands.py -i $_[0] -prefix $_[1]compara -qsub_memory 13g -q long\n", "Commands in $_[0] successfully submited!\n";
    return 1;
};

sub makeGFF
## Generates .gff file for any file with appropriate columns: chr start end
## Arguments: path_to_peak path_to_output
{
    print STDERR "Generating $_[1]\n";
    system ("awk '\$1!~\"_\"&&\$1!~\"#\"{z=\$1\":\"\$2\"-\"\$3; gsub(\"M\", \"MT\"); gsub(\"chr\", \"\"); print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"z}' $_[0] | sort -k1,1 -k2,2n | uniq > $_[1]") 
	== 0 or die ("failed to generate gff for $_[1]");
    return 1;
};

sub uniq 
## Returns non-redundant list of array elements
{
    return keys %{{ map { $_ => 1 } @_ }};
};

sub renameSpecies
## subroutine for renaming species common names -> db names for compara
## accepts a single array of species common names to be matched against 
## the species_name hash
{
    for (my $i = 0; $i <= $#_; $i++)
    {
	my @db_name_list;
	my $db_name;
	my $isFound = 0;
	for my $key (keys %species_names)
	    {
		if (($key =~ "@_[$i]" || "@_[$i]" =~ $species_names{$key})&& $isFound == 0)
		{
		    $db_name = $species_names{$key};
		    $isFound ++;
		}
		elsif (($key =~ "@_[$i]" || "@_[$i]" =~ $species_names{$key}) && $isFound > 0)
		{
		    print STDERR "Warning! more than one instance of @_[$i] found!\n";
		    print STDERR "Used organism: $db_name\t Redundant organism: $species_names{$key}\n";
		};
	    };
	if ($isFound > 0)
	{
	    @_[$i] = $db_name;
	}
	else
	{
	    die ("Couldn't find @_[$i]! Invalid Organism?")
	};
    };
    return 1;
};
