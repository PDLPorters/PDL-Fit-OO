#!perl

##########################################################################
##### checks argument parsing and error messages for sub Householder #####
##########################################################################

# If one of these tests fails, chances are really good that the
# diagnostic help in the pod is also out-of-date.


use Test::More tests => 58;
use Test::Exception;
use Test::Warn;
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

throws_ok { $coefs = $data->Householder(2) } qr/scalar as last defined argument/,
	'needs more than one argument';
like($@, qr/argument 1/, 'Error was with first argument');

throws_ok { $coefs = $data->Householder($t, 2) } qr/scalar as last defined argument/,
	'polynomial order cannot be the last argument';
like($@, qr/argument 2/, 'Error was with second argument');

throws_ok { $coefs = $data->Householder(2, 1) } qr/scalar as last defined argument/,
	'scalar is not a valid dependent piddle';
like($@, qr/argument 2/, 'Error was with second argument');

throws_ok { $coefs = $data->Householder(2, $bad_t ) } qr/incompatible dimensions/,
	'all piddles must have compatible dimensions';
like($@, qr/argument 2/, 'Error was with second argument');

throws_ok { $coefs = $data->Householder(2, 1, $t) } qr/multiple scalars/,
	'only one scalar (polynomial order) allowed';
like($@, qr/argument 2/, 'Error was with second argument');

throws_ok { $coefs = $data->Householder(2, \&PDL::sin, 1, $t) } qr/multiple scalars/,
	'only one scalar (polynomial order) allowed';
like($@, qr/argument 3/, 'Error was with third argument');

$_ = 5;
throws_ok { $coefs = $data->Householder(2, sin, $t) } qr/multiple scalars/,
	'inadvertently called Perl\'s sine function on $_';
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

throws_ok { $coefs = $data->Householder(\&my_sin) } qr/function reference as last defined argument/,
	'needs more than one argument';
like($@, qr/argument 1/, 'Error was with first argument');

throws_ok { $coefs = $data->Householder($t, \&PDL::sin) } qr/function reference as last defined argument/,
	'function reference cannot be the last argument';
like($@, qr/argument 2/, 'Error was with second argument');

throws_ok { $coefs = $data->Householder(\&PDL::sin, 1) } qr/scalar as last defined argument/,
	'scalar is not a valid dependent piddle';
like($@, qr/argument 2/, 'Error was with second argument');

throws_ok { $coefs = $data->Householder(\&PDL::sin, \&PDL::cos, $bad_t ) }
	qr/incompatible dimensions/, 'all piddles must have compatible dimensions';
like($@, qr/argument 3/, 'Error was with third argument');

my $dying_func = sub { die "I'm outa here!" };
throws_ok { $coefs = $data->Householder(\&PDL::sin, $dying_func, $t) }
	qr/Function reference must be able to process piddles/,
	'must not die when processing piddles';
like($@, qr/argument 2/, 'Error was with second argument');
throws_ok { $coefs = $data->Householder(\&sin, $t) }
	qr/Function reference must be able to process piddles/,
	'must qualify sine, cosine, etc using PDL::...';
like($@, qr/argument 1/, 'Error was with first argument');


my $weird_func1 = sub { return $bad_t };
my $weird_func2 = sub { return 5 };
throws_ok { $coefs = $data->Householder(\&PDL::sin, \&PDL::cos, $weird_func1, $t) }
	qr/Function reference must return a piddle with dimensions compatible/
	, 'function must return compatible dimensions';
like($@, qr/argument 3/, 'error was with third argument');

throws_ok { $coefs = $data->Householder(\&PDL::sin, $weird_func2, \&PDL::cos, $t) }
	qr/Function reference must return a piddle with dimensions compatible/
	, 'function must return a piddle';
like($@, qr/argument 2/, 'error was with second argument');


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
	qr/Expected a piddle but got something else/
	, 'an anonymous array is not a valid argument';
like($@, qr/argument 4/, 'Error was with the fourth argument');

throws_ok { $coefs = $data->Householder($t**2, [$t], \&PDL::sin, $t) }
	qr/Must be a scalar, function reference, or piddle/
	, 'an anonymous array is not a valid argument';
like($@, qr/argument 2/, 'Error was with the second argument');

throws_ok { $coefs = $data->Householder({}, $t**2, \&PDL::sin, $t) }
	qr/Must be a scalar, function reference, or piddle/
	, 'an anonymous hash is not a valid argument';
like($@, qr/argument 1/, 'Error was with the first argument');


# Check for undef processing and argument index reporting
# Turn off printed warnings for now
note('Check how undef effects reporting of indices of bad arguments');
PDL::Fit::Householder::no_warn_on_undef();

# First a null test to get our bearings; this test already passed above.
# All of the following attempts will croak with an incompatible
# dimensions error; our concern is what it tells us about the index of
# the argument that had the incompatible dimension:
throws_ok { $coefs = $data->Householder($bad_t, $t) }
	qr/incompatible dimensions/, 'must supply compatible piddles';
like($@, qr/argument 1/, 'Error was with the first argument');

throws_ok { $coefs = $data->Householder(undef, $bad_t, $t) }
	qr/incompatible dimensions/, 'must supply compatible piddles';
like($@, qr/argument 2/, 'Error was with the second argument');

throws_ok { $coefs = $data->Householder($bad_t, $t, undef) }
	qr/incompatible dimensions/, 'must supply compatible piddles';
like($@, qr/argument 1/, 'Error was with the first argument');

throws_ok { $coefs = $data->Householder($t, undef, $bad_t, undef, $t, undef) }
	qr/incompatible dimensions/, 'must supply compatible piddles';
like($@, qr/argument 3/, 'Error was with the third argument');

PDL::Fit::Householder::warn_on_undef();

# Check for undef processing an warning
note('Check argument indices of undef warnings');

# We'll use this wrapper function to ignore errors - we're only
# interested in warnings for this part.
sub dont_die_on_error {
	eval { $coefs = $data->Householder(@_) };
	$@ = undef;
}

warning_like { dont_die_on_error($t, undef, $bad_t) }
	{carped => qr/Argument 2/}, 'Properly reported second arg as undefined';
warnings_like { dont_die_on_error($t, undef, $bad_t, undef) }
	[{carped => qr/Argument 2/}, {carped => qr/Argument 4/}]
	, 'Properly reported second and fourth args as undefined';
my $undefined_var;
warnings_like { dont_die_on_error(undef, $undefined_var, $t, undef, $bad_t, undef) }
	[{carped => qr/Argument 1/}, {carped => qr/Argument 2/}
		,{carped => qr/Argument 4/}, {carped => qr/Argument 6/}]
	, 'Properly reported many args as undefined';

#done_testing();
