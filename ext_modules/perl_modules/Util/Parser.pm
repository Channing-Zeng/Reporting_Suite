package Util::Parser;

=head2 column_hash

  Usage   : $hash = Util::Parser::column_hash(-file => ..., -kcols => ..., -vcols => ..., ...)
  Function: Parse a file and construct a hash using given columns
  Args    : Named Parameter:
  	    -file	The input file
	    -kcols	The columns (start from 1) for key.  Multiple columns are 
	    	        separated by :.  Default to 1
	    -vcols	The columns (start from 1)for value.  Multiple columns are 
	    		separated by :.  Default to 2
	    -dem	The delimiter to get columns.  Default to "\t"
	    -jdem	The delimiter to join multiple columns for key and value.
	    		If not specified, will use the "-dem".

=cut

sub column_hash {
    my %pm = @_;
    my @kcols = $pm{ "-kcols" } ? split( /:/, $pm{ "-kcols" } ) : (1);
    my @vcols = $pm{ "-vcols" } ? split( /:/, $pm{ "-vcols" } ) : (2);
    @kcols = map { $_ - 1; } @kcols;
    @vcols = map { $_ - 1; } @vcols;
    my $dem = $pm{ "-dem" } ? $pm{ "-dem" } : "\t";
    my $jdem = $pm{ "-jdem" } ? $pm{ "-jdem" } : $dem;

    open( IN, $pm{ "-file" } ) || die "Can't open file $pm{ -file }";
    my %hash;
    while( <IN> ) {
        chomp;
	my @a = split( $dem, $_, -1 );
	my $k = join( $jdem, @a[@kcols] );
	my $v = join( $jdem, @a[@vcols] );
	$hash{ $k } = $v;
    }
    close( IN );
    return \%hash;
}

=head2 range

  Usage   : @index = range($str);
  Function: Convert a range specifier (1-based) to an array of index (0-based).
 	    Multiple, non-consecutive ranges are concatenated with ":".  Three types of
	    ranges are accepted:
	    1. A single number.  Use as is.
	    2. Consecutive range can be specified using .. operator.  For example 1..10
	    3. Evenly spaced range.  For example, 1..10s2 will return (0, 2, 4, 6, 8).

	    Program dies if the range is illegal.

  Argument: A string as described above
  Returns : An index array (0 - based).

=cut

sub range {
    my $c = shift;
    my @c = ();
    my @a = split( /:/, $c );
    foreach(@a) {
        if ( /^\d+$/ ) {
            push( @c, $_ );
        } elsif ( /^(\d+)\.\.(\d+)s(\d+)$/ ) {
            for(my $i = $1; $i <= $2; $i += $3) {
                push( @c, $i );
            }
        } elsif ( /^(\d+)\.\.(\d+)$/ ) {
            push( @c, ($1 .. $2) );
        } else {
            die "Illegal range '$_'.  Please use 'perldoc Util::Parser' for documents.";
        }
    }
    @c = map { $_-1; } @c;
    return @c;
}

=head2 indexes

  Usage   : @index = indexes($str);
  Function: Convert a indexes specifier (1-based) to an array of indexes (0-based).
 	    Multiple, non-consecutive ranges are concatenated with ":".  Three types of
	    ranges are accepted:
	    1. A single number.  Use as is.
	    2. Consecutive range can be specified using .. operator.  For example 1..10
	    3. Evenly spaced range.  For example, 1..10s2 will return (0, 2, 4, 6, 8).

	    Program dies if the range is illegal.

  Argument: A string as described above
  Returns : An index array (0 - based).

=cut

sub indexes {
    my $c = shift;
    my @c = ();
    my @a = split( /:/, $c );
    foreach(@a) {
        if ( /^\d+$/ ) {
            push( @c, $_ );
        } elsif ( /^(\d+)\.\.(\d+)s(\d+)$/ ) {
            for(my $i = $1; $i <= $2; $i += $3) {
                push( @c, $i );
            }
        } elsif ( /^(\d+)\.\.(\d+)$/ ) {
            push( @c, ($1 .. $2) );
        } else {
            die "Illegal range '$_'.  Please use 'perldoc Util::Parser' for documents.";
        }
    }
    @c = map { $_-1; } @c;
    return @c;
}

return 1;
