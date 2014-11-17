package Util::Nexus;
use strict;

=head1 NAME

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 EXAMPLES

=head1 AUTHOR

  Zhongwu Lai

=cut

=head1 METHODS

=cut

use strict;

=head2 new

  Usage   :
  Function:
  Args    :
  Returns : An object reference

=cut

sub new {
    my $proto = shift;
    my $class = ref( $proto ) || $proto;

    my $self = {};

    my %p = @_;
    while( my ($k, $v) = each %p ) {
        $k = lc( $k );
        $k =~ s/^-/_/;
        $self->{ $k } = $v;
    }
    bless $self, $class;
    return $self;
}

=head1 COPYRIGHT

  Copyright (c) 2007, AstraZeneca.  All Rights Reserved.

=cut

1;
