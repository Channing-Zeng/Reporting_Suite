#!/usr/bin/perl -w

use lib '/users/kdld047/lib/perl5';
use Getopt::Std;
use Stat::Basic;
use strict;

our ($opt_i, $opt_c, $opt_v, $opt_h, $opt_H, $opt_I, $opt_f, $opt_V, $opt_S);

Usage() unless(getopts( 'hHSi:c:v:I:f:V:' ));
Usage() if ( $opt_H );
my $stat = new Stat::Basic;

$opt_i = $opt_i ? $opt_i : 1;
$opt_c = $opt_c ? $opt_c : 2;
$opt_v = $opt_v ? $opt_v : 3;

my @i = getcols( $opt_i );
my @c = getcols( $opt_c );
my @v = getcols( $opt_v );
my $V = $opt_V;

@i = map { $_ - 1; } @i;
@c = map { $_ - 1; } @c;
@v = map { $_ - 1; } @v;

my $id = join("\t", (map { ""; } @i));
my $hdr;
my @hdrs;
if ( $opt_h ) {
    $hdr = <>; chomp( $hdr );
    $hdr =~ s/\r//g;
    @hdrs = split( /\t/, $hdr );
    $id = join("\t", @hdrs[@i]);
}
if ( $opt_I ){
    $id = $opt_I;
    $id =~ s/\\t/\t/g;
}
my %c;  # for categories

my %hash;
while( <> ) {
    chomp;
    s/\r//g;
    my @a = split( /\t/ );
    my $k = join("\t", @a[@i]);
    my $c = join(",", @a[@c]);
    if ( $opt_V ) {
        my $ck = $opt_h ? $c : $c;
        push( @{ $hash{ $k }->{ $ck } }, $V );
	$c{ $ck } = 1;
    } else {
	foreach(@v) {
	    my $v = $a[$_];
	    #my $ck = $opt_h ? "$c-$hdrs[$_]" : "$c";
	    my $ck = $opt_h ? "$c" : "$c";
	    push( @{ $hash{ $k }->{ $ck } }, $v) if ( defined($v) );
	    $c{ $ck } = 1;
	}
    }
}

my @cats = sort(keys %c);

if ( @v > 1 ) {
    my @tmp;
    if ( $opt_h ) {
        @tmp = @hdrs[@v];
    } else {
        @tmp = map { $_ + 1; } @v;
    }
    my @thdr = ($id);
    foreach my $c (@cats) {
        foreach my $t (@tmp) {
	    push(@thdr, "$c-$t");
	}
    }
    print join("\t", @thdr), "\n";
} else {
    print join("\t", $id, @cats), "\n";
}
if ( $opt_S ) {
    my @ids = sort(keys %hash);
    foreach my $k (@ids) {
	my $v = $hash{ $k };
	my @t = map { $v->{ $_ } ? (($opt_f && @{ $v->{ $_ } }>1) ? $stat->$opt_f( $v->{ $_ } ) : (@v > 1 ? join("\t", @{ $v->{ $_ } }) : join( ";", @{ $v->{ $_ } }) )) : ""; } @cats;
	print join("\t", $k, @t), "\n";
    }
} else {
    while( my ($k, $v) = each %hash ) {
	#print STDERR "K: $k\n";
	my @t = map { $v->{ $_ } ? (($opt_f && @{ $v->{ $_ } }>1) ? $stat->$opt_f( $v->{ $_ } ) : (@v > 1 ? join("\t", @{ $v->{ $_ } }) : join( ";", @{ $v->{ $_ } }) )) : ""; } @cats;
	print join("\t", $k, @t), "\n";
    }
}

sub getcols {
    my $c = shift;
    my @c = ();
    my @a = split( /:/, $c );
    foreach(@a) {
        if ( /^\d$/ ) {
            push( @c, $_ );
        } elsif ( /^(\d+)\.\.(\d+)s(\d+)$/ ) {
            for(my $i = $1; $i <= $2; $i += $3) {
                push( @c, $i );
            }
        } elsif ( /^(\d+)\.\.(\d+)$/ ) {
            push( @c, ($1 .. $2) );
        } else {
            push( @c, $_ );
        }
    }
    return @c;
}

sub Usage {

    print STDERR<<USAGE;

$0 [-hHS] [-i ID_columns] [-c category_columns] [-v value_columns] [-I new_id_column_name] [-f function] [-V VALUE] file

The program will pivot a tall table into a fat table.

The following options are available:

-i ID_columns for row in the matrix.  Column starts from 1.  Can use multiple columns, separated by : 
   or using .. combined with 's', such as 1..10s2.  Default to 1st column.
-c Category columns for column in the matrix.  Default to 2nd column.
-v Value columns.  Default to 3rd column.
-V The "VALUE", such as "Yes", "Mut" or "1" to indicate the key is observed.  Should be used rarely.
-I The name for new ID column if it's different.
-h Indicate whether header exists.  Default no.
-f The function to calculate the value when there're duplicates, such as mean, max, min, abmax, and abmin.
-S Sort the ID column
-H Usage
USAGE
    exit();
}
