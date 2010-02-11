#!perl

use Test::More;
use Test::Exception;
#use Test::More tests => 23;
use PDL;
use PDL::Fit::Householder;

use strict;
use warnings;

# Build a little bit of data
my $t = sequence(100) / 10;
my $bad_t = sequence(50) / 10;
my ($a, $b) = (4, -2);
my $data = $a * sin($t) + $b * cos($t) + grandom($t) * 0.1;
my $coefs;

##### check argument parsing for sub Householder #####

# Check basic requirements of incoming arguments:
note('Basic requirements');
throws_ok { $coefs = $data->Householder } qr/Householder\(\)/,
	'requires at least two arguments';
throws_ok { $data->Householder } qr/Householder\(\)/,
	'does not allow void context';
throws_ok { PDL::Fit::Householder::Householder('string', $t) } qr/Householder\(\)/,
	'must act upon a piddle';


# Check polynomial processing
note('Polynomial errors');

throws_ok { $coefs = $data->Householder(2) } qr/scalar as last argument/,
	'needs more than one argument';
like($@, qr/last argument/, 'Error was with last argument');

throws_ok { $coefs = $data->Householder($t, 2) } qr/scalar as last argument/,
	'polynomial order cannot be the last argument';
like($@, qr/last argument/, 'Error was with last argument');

throws_ok { $coefs = $data->Householder(2, 1) } qr/scalar as last argument/,
	'scalar is not a valid dependent piddle';
like($@, qr/last argument/, 'Error was with last argument');

throws_ok { $coefs = $data->Householder(2, $bad_t ) } qr/incompatible dimensions/,
	'all piddles must have compatible dimensions';
like($@, qr/last argument/, 'Error was with last argument');

throws_ok { $coefs = $data->Householder(2, 1, $t) } qr/multiple scalars/,
	'only one scalar (polynomial order) allowed';
like($@, qr/argument 2/, 'Error was with second argument');

throws_ok { $coefs = $data->Householder(2.2, $t) } qr/integer/,
	'polynomial order must be an integer';
like($@, qr/argument 1/, 'Error was with first argument');

throws_ok { $coefs = $data->Householder(-1, $t) } qr/negative/,
	'polynomial order cannot be negative';
like($@, qr/argument 1/, 'Error was with first argument');


# Check function processing
note('Anonymous functions and function reference errors');

sub my_sin { return $_[0]->sin() }

throws_ok { $coefs = $data->Householder(\&my_sin) } qr/function reference as last argument/,
	'needs more than one argument';
like($@, qr/last argument/, 'Error was with last argument');

throws_ok { $coefs = $data->Householder($t, \&PDL::sin) } qr/function reference as last argument/,
	'function reference cannot be the last argument';
like($@, qr/last argument/, 'Error was with last argument');

throws_ok { $coefs = $data->Householder(\&PDL::sin, 1) } qr/scalar as last argument/,
	'scalar is not a valid dependent piddle';
like($@, qr/last argument/, 'Error was with last argument');

throws_ok { $coefs = $data->Householder(\&PDL::sin, \&PDL::cos, $bad_t ) }
	qr/incompatible dimensions/, 'all piddles must have compatible dimensions';
like($@, qr/last argument/, 'Error was with last argument');

my $dying_func = sub { die "I'm outa here!" };
throws_ok { $coefs = $data->Householder(\&PDL::sin, $dying_func, $t) }
	qr/Function reference must be able to process piddles/,
	'must not die when processing piddles';
like($@, qr/argument 2/, 'Error was with second argument');

my $weird_func = sub { return $bad_t };
throws_ok { $coefs = $data->Householder(\&PDL::sin, \&PDL::cos, $weird_func, $t) }
	qr/return .* compatible/, 'must return compatible dimensions';
like($@, qr/argument 3/, 'error was with third argument');


# Check piddle processing
note('Arbitrary piddle arguments and their errors');

throws_ok { $coefs = $data->Householder(sin($bad_t), $t, $bad_t) }
	qr/incompatible dimensions/, 'must supply compatible piddles';
like($@, qr/argument 1/, 'Error was with the first argument');

throws_ok { $coefs = $data->Householder(sin($t), $t, $bad_t) }
	qr/incompatible dimensions/, 'must supply compatible piddles, round 2';
like($@, qr/argument 3/, 'Error was with the third argument');


# Check other errors
note('Other input errors');
throws_ok { $coefs = $data->Householder($t**2, \&PDL::sin, $t, [$t]) }
	qr/Expected a piddle but got something else/, 'an anonymous array is not a valid argument';
like($@, qr/argument 2/, 'Error was with the second argument');

throws_ok { $coefs = $data->Householder($t**2, [$t], \&PDL::sin, $t) }
	qr/Must be/, 'an anonymous array is not a valid argument';
like($@, qr/argument 2/, 'Error was with the second argument');

throws_ok { $coefs = $data->Householder({}, $t**2, \&PDL::sin, $t) }
	qr/Must be/, 'an anonymous hash is not a valid argument';
like($@, qr/argument 1/, 'Error was with the first argument');


# Check that pdl(1) works in creating a piddle of constants
# Check undef handling with and without warnings enabled

# 2) check _check_piddle to make sure it works right
# 3) check that _backsub works correctly by constructing a matrix and y,
#    finding coefs, and then making sure that $matrix x $coefs == $y
# 4) check that _householder gives the same results as PDL::MatrixOps
#    for diagonal matrices
# 5) Check that generated signals with noise can be fit with reasonable
#    results.

done_testing();
