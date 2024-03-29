=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=head1 NAME

Bio::EnsEMBL::Hive::RunnableDB::DumpMultiAlign::Compress.pm

=head1 SYNOPSIS

This RunnableDB module is part of the DumpMultiAlign pipeline.

=head1 DESCRIPTION

This RunnableDB module runs gzip -9

=cut


package Bio::EnsEMBL::Compara::RunnableDB::DumpMultiAlign::Compress;

use strict;
use base ('Bio::EnsEMBL::Compara::RunnableDB::BaseRunnable');

sub fetch_input {
    my $self = shift;

}

sub run {
    my $self = shift;

    #
    #Run gzip -9 command (with force option)
    #
    my $output_file = $self->param('output_dir') . "/" . $self->param('output_file');

    #Check existence of output_file
    return unless (-e $output_file);

    my $cmd = "gzip -f -9 " . $output_file;
    if(my $return_value = system($cmd)) {
        $return_value >>= 8;
        die "system( $cmd ) failed: $return_value";
    }

    #
    #If maf_output_dir defined, move maf file from emf directory to maf
    #directory
    #
    if ($self->param('maf_output_dir')) {
	my $mv_cmd = "mv " . $self->param('output_dir') . "/" . $self->param('output_file') . ".gz " . $self->param('maf_output_dir');
	if(my $return_value = system($mv_cmd)) {
	    $return_value >>= 8;
	    die "system( $cmd ) failed: $return_value";
	}
    }

}

sub write_output {
    my $self = shift @_;

}

1;
