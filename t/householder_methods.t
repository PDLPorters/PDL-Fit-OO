######################################
##### Check numerical algorithms #####
######################################

use Test::More;
use Test::LectroTest;
use PDL;
use PDL::Fit::Householder;
use PDL::NiceSlice;

use strict;
use warnings;

sub backsub_test {
	my ($matrix_diags, $coefs) = @_;

	# build an upper-triangular matrix with the given diagonal and an
	# otherwise random upper triangle:
	my $matrix = stretcher($matrix_diags);
	my $upper_triangle = $matrix->where($matrix->xvals > $matrix->yvals);
	$upper_triangle .= random($upper_triangle);

	# Compute y, then use _backsub to solve for $coefs
	my $y = $matrix x $coefs->dummy(0);
	my $backsub_coefs = _backsub($y, $matrix)->squeeze;
	return all(approx($backsub_coefs, $coefs, 1e-8));
}


# Test 1
Property {
	##[
		a <- Float(range=>[-100, 100]),
		b <- Float(range=>[-100, 100]),
		c <- Float(range=>[-100, 100]),
		d <- Float(range=>[-100, 100]),
	]##

	# Make sure our matrix is not singular (to machine precision):
	return $tcon->retry if ($a + $b == $a or $a + $b == $b);

	backsub_test(pdl($a, $b), pdl($c, $d));
}, name => '_backsub works for 2x2 matrices';

# Test 2
Property {
	##[
		a <- Float(range=>[-100, 100]),
		b <- Float(range=>[-100, 100]),
		c <- Float(range=>[-100, 100]),
		d <- Float(range=>[-100, 100]),
		e <- Float(range=>[-100, 100]),
		f <- Float(range=>[-100, 100]),
	]##

	# Make sure our matrix is not singular:
	return $tcon->retry if ($a == 0 or $b == 0 or $c == 0);

	backsub_test(pdl($a, $b, $c), pdl($d, $e, $f));
}, name => '_backsub works for 3x3 matrices';


sub householder_test {
	my ($matrix_diags, $coefs) = @_;

	# build an upper-triangular matrix with the given diagonal and an
	# otherwise random upper triangle:
	my $matrix = stretcher($matrix_diags);
	my $upper_triangle = $matrix->where($matrix->xvals > $matrix->yvals);
	my $lower_triangle = $matrix->where($matrix->xvals < $matrix->yvals);
	$lower_triangle .= $upper_triangle .= random($upper_triangle);

	# Compute y, then use _backsub to solve for $coefs
	my $y = $matrix x $coefs->dummy(0);
	$y->_Householder($matrix);
	my $Householder_coefs = _backsub($y, $matrix)->squeeze;
	return all(approx($Householder_coefs, $coefs, 1e-8));
}

# Test 3
Property {
	##[
		a <- Float(range=>[-100, 100]),
		b <- Float(range=>[-100, 100]),
		c <- Float(range=>[-100, 100]),
		d <- Float(range=>[-100, 100]),
	]##

	# Make sure our matrix is not singular (to machine precision):
	return $tcon->retry if ($a + $b == $a or $a + $b == $b);

	householder_test(pdl($a, $b), pdl($c, $d));
}, name => '_Householder works for symmetric 2x2 matrices';

# Test 4
Property {
	##[
		a <- Float(range=>[-100, 100]),
		b <- Float(range=>[-100, 100]),
		c <- Float(range=>[-100, 100]),
		d <- Float(range=>[-100, 100]),
		e <- Float(range=>[-100, 100]),
		f <- Float(range=>[-100, 100]),
	]##

	# Make sure our matrix is not singular:
	return $tcon->retry if ($a == 0 or $b == 0 or $c == 0);

	householder_test(pdl($a, $b, $c), pdl($d, $e, $f));
}, name => '_Householder works for symmetric 3x3 matrices';

my $t = sequence(100) / 10;
my $eps = 0.01;
my $noise_factor = 0.001;

# Test 5
Property {
	##[
		a <- Float(range=>[-100, 100]),
		b <- Float(range=>[-100, 100]),
	]##

	my $signal = $a + $b * $t + grandom($t) * $noise_factor;
	my $given_coefs = pdl($a, $b);
	my $coefs = $signal->Householder(1, $t);
	$tcon->note("Deltas are " . join('; ', ($coefs - $given_coefs)->list));
	return all(approx($coefs, $given_coefs, $eps));
}, name => 'Linear (first order polynomial) fitting works';

# Test 6
Property {
	##[
		a <- Float(range=>[-100, 100]),
		b <- Float(range=>[-100, 100]),
		c <- Float(range=>[-100, 100]),
	]##

	my $signal = $a + $b * $t + $c * $t**2 + grandom($t) * $noise_factor;
	my $given_coefs = pdl($a, $b, $c);
	my $coefs = $signal->Householder(2, $t);
	$tcon->note("Deltas are " . join('; ', ($coefs - $given_coefs)->list));
	return all(approx($coefs, $given_coefs));
}, name => 'Second order polynomial fitting works';

# Test 7
Property {
	##[
		a <- Float(range=>[-100, 100]),
		b <- Float(range=>[-100, 100]),
	]##

	my $signal = $a * sin($t) + $b * cos($t) + grandom($t) * $noise_factor;
	my $given_coefs = pdl($a, $b);
	my $coefs = $signal->Householder(\&PDL::sin, \&PDL::cos, $t);
	$tcon->note("Deltas are " . join('; ', ($coefs - $given_coefs)->list));
	return all(approx($coefs, $given_coefs));
}, name => 'Supplied function fitting works';

# Test 8
Property {
	##[
		a <- Float(range=>[-100, 100]),
		b <- Float(range=>[-100, 100]),
		c <- Float(range=>[-100, 100]),
	]##

	my $signal = $a * sin($t) + $b + $c * $t + grandom($t) * $noise_factor;
	my $given_coefs = pdl($a, $b, $c);
	my $coefs = $signal->Householder(\&PDL::sin, 1, $t);
	$tcon->note("Deltas are " . join('; ', ($coefs - $given_coefs)->list));
	return all(approx($coefs, $given_coefs));
}, name => 'Combined polynomial and supplied function fitting works';

# Test 9
Property {
	##[
		a <- Float(range=>[-100, 100]),
		b <- Float(range=>[-100, 100]),
		c <- Float(range=>[-100, 100]),
	]##

	my $signal = $a * sin($t) + $b + $c * $t + grandom($t) * $noise_factor;
	my $sine = sin($t);
	my $given_coefs = pdl($a, $b, $c);
	my $coefs = $signal->Householder($sine, 1, $t);
	$tcon->note("Deltas are " . join('; ', ($coefs - $given_coefs)->list));
	return all(approx($coefs, $given_coefs));
}, name => 'Combined piddle and polynomial fitting works';

# Test 10
Property {
	##[
		a <- Float(range=>[-100, 100]),
		b <- Float(range=>[-100, 100]),
		c <- Float(range=>[-100, 100]),
		d <- Float(range=>[-100, 100]),
	]##

	my $signal = $a * sin($t) + $b + $c * $t + $d * exp($t) + grandom($t) * $noise_factor;
	my $sine = sin($t);
	my $given_coefs = pdl($a, $b, $c, $d);
	my $coefs = $signal->Householder($sine, 1, \&PDL::exp, $t);
	$tcon->note("Deltas are " . join('; ', ($coefs - $given_coefs)->list));
	return all(approx($coefs, $given_coefs));
}, name => 'Everything combined works';



# 5) Check that generated signals with noise can be fit with reasonable
#    results.

#done_testing();
