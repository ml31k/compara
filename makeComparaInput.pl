#!/usr/bin/perl

## Michael - 2013.01.29

## This script is designed to function as a wrapper for the compara package.
## It will build the directory structure and appropriately name files for
## subsequent analysis.
## It will generate text files with commands that can subsequently be submitted
## to the cluster

use strict;
use Getopt::Long;
use Cwd;

my $cwd = getcwd();
my $inputFile;
my $help;
my $path = "/home/mliang/mdwilson/Compara_2014_01_31/scripts";
my $out_species;
my %species_names;
my $is_gff;
my $alignment = "EPO";
my $mmr = "0:1";
my $do_paralogs = "FALSE";
my $auto_queue;
my $command_prefix = " ";
my $out_species_string;

GetOptions
    ("h|help" => \$help,
     "i=s" => \$inputFile,
     "o=s" => \$out_species,
     "is_gff" => \$is_gff,
     "a=s" => \$alignment,
     "mmr=s" =>\$mmr,
     "p=s" => \$do_paralogs,
     "q=s" => \$auto_queue,
     "cmd=s" => \$command_prefix
     )
    or die ("invalid commandline args\n");

my $helpmsg = join ("\n", 
		    "##############################################",
		    "Compara wrapper version 1.00 by Mimi31k",
		    "##############################################",
		    "Usage: -i inputfile [OPTIONS]",
		    "",
		    "Mandatory Arguments:",
		    "\t-i input\t\tA tab delimited file with columns Species ID(ex. TF name) path/to/peak",
		    "\t\t\t\tSpecies names must correspond to species_db.txt found in the compara directory",
		    "\t\t\t\tPeak files may be in any format with first 3 columns chr# Start End (no header)",
		    "\t\t\t\tNon .gff files will be converted and stored locally. gff files will be linked.",
		    "",
		    "Optional Arguments:",
		    "\t-h | -help\t\tOpen this help menu and exit",
		    "\t-o out_species\t\tpath/to/file containing output species. Equal to input species by default",
		    "\t-is_gff\t\t\tPeak files are in .gff format and do not require conversion. False by default",
		    "\t-a alignment\t\tAlignment mode to run Compara [EPO|PECAN]. EPO by default",
		    "\t-mmr max:min\t\tMax:Min ratios for peak overlaps. Default 0:1",
		    "\t-p do_paralogs\t\t[TRUE] Determines if R component runs paralog analysis. FALSE by default",
		    "\t-q auto_queue\t\t[FALSE] Determines if compara perl commands are sent to qsub. TRUE by default",
		    "\t-cmd command prefix\tCommand prefix for all compara commands (essential source commands, etc.)\n"		    
    );

if ($help || (not defined $inputFile) || ($alignment && $alignment ne "EPO" && $alignment ne "PECAN") || ($do_paralogs && $do_paralogs ne "TRUE" &&$do_paralogs ne "FALSE") || ($auto_queue && $auto_queue ne "FALSE"))
{
    print $helpmsg;
    exit;
};
## Auto queueing of commands is set on by default, unless specified as "FALSE"
$auto_queue = $command_prefix unless $auto_queue == "FALSE";

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
    open OUTPUT_SPECIES, $out_species or die ("couldn't open out_species file $out_species!");
    my @output_species;
    while (<OUTPUT_SPECIES>)
    {
	chomp();
	push (@output_species,$_.";");
    };
    close (OUTPUT_SPECIES);
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
	push (@{$hash_by_IDs{$a_ID}}, [$a_species, $a_peakFile]);
    }
    else 
    {
	$hash_by_IDs{$a_ID} = [];
	push (@{$hash_by_IDs{$a_ID}}, [$a_species, $a_peakFile]);
    };
};
close (INPUT);

## Setup subdirectory structure and generates peaks if necessary.
## Also prints the appropriate commands to run the perl and R portions of compara
for my $ID (keys %hash_by_IDs)
{
    if (-d $ID)
    {
	print STDERR "directory $ID exists!\n";
    }
    else
    {
	system ("mkdir $ID") == 0 or die ("failed to make subdirectory for $ID\n");
	system ("mkdir $ID/peaks") == 0 or die ("failed to make subdirectory for $ID\n");
    };
    open CMD_OUT, ">$ID/run_compara_cmd.txt" or die("failed to write CMD_OUT commands\n");
    open R_INPUT, ">$ID/rscript_input.txt" or die("failed to write R_INPUT commands\n");
    open R_INSTRUCT, ">$ID/rscript_instructions.txt" or die("failed R_INSTRUCT to write commands\n");
    open R_CMD, ">$ID/rscript_run.R" or die("failed to write R_CMD commands\n");
    for my $ID_list (\@{$hash_by_IDs{$ID}})
    {
	for my $spec_and_file (@{$ID_list}) 
	{
	    my $peakFileIn = $spec_and_file -> [1];
	    my $fileOut =  $ID."_".$spec_and_file->[0];      
	    
	    ## Link peaks if gff, else generate gff files. named: Species_db_name_ID.gff
	    if ($is_gff)
	    {
		system ("ln -s $peakFileIn $ID/peaks/$fileOut.gff") == 0 or die ("failed to link $peakFileIn\n"); 
	    }
	    else
	    {
		&makeGFF ($peakFileIn, "$ID/peaks/$fileOut.gff");
	    };
	    
	    ## Generate command to run compara:
	    print CMD_OUT ("$command_prefix ",
			   "cd $cwd ; $path/compara_match_peaks.pl ",
			   "-i $cwd/$ID/peaks/$fileOut.gff ",
			   "-o $cwd/$ID/$fileOut.compara.regions ",
			   "-s $spec_and_file->[0] ",
			   "-r $out_species_string ",
			   "-a $alignment\n");
	    ## Generate the input file list for R component
	    print R_INPUT "$spec_and_file->[0]\t$cwd/$ID/$fileOut.compara.regions\t$cwd/$ID/peaks/$fileOut.gff\n";
	};
	## Generate the instruction file for R script
	print R_INSTRUCT "$cwd/$ID/rscript_input.txt\t$ID\t$mmr\n";
	## Generate the actual R script to be executed
	print R_CMD ("source (\"$path/peakConverter.R\")\n",
		     "source (\"$path/process_compara.R\")\n",
		     "source (\"$path/batchCompara.R\")\n",
		     "do_paralogs <- $do_paralogs\n",
		     "batchComparaData(\"$cwd/$ID/rscript_instructions.txt\")\n");
	
	close (R_INPUT);
	close (R_INSTRUCT);
	close (R_CMD);
	close (CMD_OUT);
	## Submits perl commands to the cluster unless auto_queue was specified as "FALSE"
	&queue_commands("$ID/run_compara_cmd.txt", $ID) unless $auto_queue;    
    };
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
    system ("queue_commands.py -i $_[0] -prefix $_[1]compara -qsub_memory 13g") == 0 or die ("Error in submitting jobs to cluster!\n");
    print "System command: queue_commands.py -i $_[0] -prefix $_[1]compara -qsub_memory 13g\n", "Commands in $_[0] successfully submited!\n";
    return 1;
};

sub makeGFF
## Generates .gff file for any file with appropriate columns: chr start end
## Arguments: path_to_peak path_to_output
{
    print STDERR "Generating $_[1]\n";
    system ("awk '\$1!~\"_\"&&\$1!~\"#\"{z=\$1\":\"\$2\"-\"\$3; gsub(\"chr\", \"\"); print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"z}' $_[0] > $_[1]") 
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
