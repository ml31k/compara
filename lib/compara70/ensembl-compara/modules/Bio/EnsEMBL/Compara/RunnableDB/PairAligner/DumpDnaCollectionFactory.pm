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

Bio::EnsEMBL::Compara::RunnableDB::PairAligner::DumpDnaCollectionFactory

=cut

=head1 SYNOPSIS

=cut

=head1 DESCRIPTION


=cut

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Compara::RunnableDB::PairAligner::DumpDnaCollectionFactory;

use strict;
use Time::HiRes qw(time gettimeofday tv_interval);
use Bio::EnsEMBL::Analysis::Runnable::Blat;
use Bio::EnsEMBL::Analysis::RunnableDB;

use base ('Bio::EnsEMBL::Compara::RunnableDB::BaseRunnable');
use File::Path;

sub fetch_input {
  my( $self) = @_;

  if ($self->param('dna_collection_name')) {
      $self->param('collection_name', $self->param('dna_collection_name'));
  }

  die("Missing dna_collection_name") unless($self->param('collection_name'));
  die("Must specifiy dump_min_size") unless ($self->param('dump_min_size'));

  return 1;
}



sub run
{
  my $self = shift;

   return 1;
}


sub write_output {
  my( $self) = @_;

 if ($self->param('dump_nib')) {
      $self->dumpNibFilesFactory;
  }
  if ($self->param('dump_dna')) {
      $self->dumpDnaFilesFactory;
  }

  return 1;
}


##########################################
#
# internal methods
#
##########################################

sub dumpNibFilesFactory {
  my $self = shift;

  my $starttime = time();
  my $dna_collection = $self->compara_dba->get_DnaCollectionAdaptor->fetch_by_set_description($self->param('collection_name'));
  my $dump_loc = $dna_collection->dump_loc;

  unless (defined $dump_loc) {
    die("dump_loc directory is not defined, can not dump nib files\n");
  }
  foreach my $dnafrag_chunk_set (@{$dna_collection->get_all_DnaFragChunkSets}) {
      my $output_id;
      foreach my $dnafrag_chunk (@{$dnafrag_chunk_set->get_all_DnaFragChunks}) {
          next if ($dnafrag_chunk->length <= $self->param('dump_min_size'));

          my $nibfile = "$dump_loc/". $dnafrag_chunk->dnafrag->name . ".nib";
          
          #don't dump nibfile if it already exists
          next if (-e $nibfile);
          
          $output_id->{'DnaFragChunk'} = $dnafrag_chunk->dbID;
          $output_id->{'collection_name'} = $self->param('collection_name');
          $output_id->{'dump_loc'} = $self->param('dump_loc');
          $output_id->{'genome_db_id'} = $self->param('genome_db_id');
          
          #Add dataflow to branch 2
          $self->dataflow_output_id($output_id,2);
          
      }
  }

  if($self->debug){printf("%1.3f secs to dump nib for \"%s\" collection\n", (time()-$starttime), $self->param('collection_name'));}

  return 1;
}

sub dumpDnaFilesFactory {
  my $self = shift;

   my $dna_collection = $self->compara_dba->get_DnaCollectionAdaptor->fetch_by_set_description($self->param('collection_name'));

  foreach my $dnafrag_chunk_set (@{$dna_collection->get_all_DnaFragChunkSets}) {
      my $type;
      my $output_id;

      #Only dump single dnafrag_chunk if it is larger than the dump_min_size - otherwise it gets stored in the database
      my $dnafrag_chunks = $dnafrag_chunk_set->get_all_DnaFragChunks;
      next if (@$dnafrag_chunks == 1 && $dnafrag_chunks->[0]->length <= $self->param('dump_min_size'));

      $type = "DnaFragChunkSet";
      $output_id->{$type} = $dnafrag_chunk_set->dbID;
      $output_id->{'collection_name'} = $self->param('collection_name');

      #Add dataflow to branch 2
      $self->dataflow_output_id($output_id,2);

  }
   return 1;
}

1;
