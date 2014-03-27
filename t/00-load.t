#!perl -T
use 5.006;
use strict;
use warnings FATAL => 'all';
use Test::More;

plan tests => 3;

BEGIN {
    use_ok( 'Stats' ) || print "Bail out!\n";
    use_ok( 'Plot' ) || print "Bail out!\n";
    use_ok( 'Forge' ) || print "Bail out!\n";
}

diag( "Testing Stats $Stats::VERSION, Perl $], $^X" );
