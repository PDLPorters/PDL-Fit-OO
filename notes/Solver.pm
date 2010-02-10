package PDL::Fit::Solver;

use warnings;
use strict;
use Carp;

use PDL;
use PDL::Core ':Internal';
use PDL::NiceSlice;
#use PDL::Iterator;

our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
use version;
$VERSION = qv('0.0.1');

use PDL::Exporter;
@ISA = qw(PDL::Exporter);
@EXPORT_OK = qw(solve Householder);
%EXPORT_TAGS = (
	Func => [@EXPORT_OK],
	Internal => [qw/Givens Householder SVD build_anon_func build_anon_poly/]
);

our $DEBUG = 1;

sub Givens {
	
}

=head2 Householder

=for ref

Performs a Householder diagonalization of the supplied matrix.

=for usage

 Usage:
   PIDDLE = PIDDLE->Householder(PIDDLE, OPTIONS)
   HASH = PIDDLE->Householder(PIDDLE, OPTIONS)

=for options

At the moment, OPTIONS is optional.  If set to a true value, it will perform
the calculations in place.  If not specified or explicitly set to a false
value, it will make a copy of the given piddle and matrix piddle.

The piddles you supply can have the following arguments:

 self   A      Comments     
 ============================================================================
 1-dim  2-dim  Treats $self like $self->dummy(0); $A->dim(0) == $self->dim(0)
 2-dim  2-dim  $A->dim(0) == $self->dim(1)

=cut

sub Householder {
	my ($y, $A, $inplace) = @_;
	
	# Basic sanity checks:
	$A->ndims > 1 or barf("Usage: PDL::Householdr");	# $A must be meaningful
	$y->ndims > 0 or barf("Usage: PDL::Housholder");	# $y must be meaningful
	$A->dim(0) > $A->dim(1) or barf("Usage: PDL::Householder");
	
	# Dimensions must agree:
	$y = $y->dummy(0) if ($y->ndims == 1);
	$y->dim(1) == $A->dim(0) or barf("Usage: PDL::Householder");
	
	my $original_y;
	# Copy data if necessary
	if ($inplace) {
		# We only need the original y values if we're calculating
		# chi-squared, which only happens when the want an array
		# returned
		$original_y = $y->copy if wantarray;
	}
	else {
		$original_y = $y;
		($y, $A) = ($y->copy, $A->copy);
	}

	$y->_Householder_simple($A);
	my $coefs = $y->_backsub($A);
	
	# Eventually return coefficients, residual, and method
	return $coefs if (defined wantarray and not wantarray);
	my ($residual);
	return (
		coefs => $coefs,
		residual => $residual,
		method => 'Householder',
		diagonalized_matrix => $A,
		manipulated_results => $y,
	) if (wantarray);
}

use Inline (Pdlpp => Config => AUTO_INCLUDE => '#include <math.h>');
no PDL::NiceSlice;

# Only use this more complicated construction if I must:
#		Pars => 'y(l, m); A(n, m); [o] coefs(l, n);',
use Inline
	Pdlpp => q{
	pp_def('_Householder_fast',
		Pars => 'y(i, m); A(n, m); [t]v(m);',
		Doc => q{
Note that this function is supposed to be just like C<_Householder_simple>,
but a little more complicated and a little bit faster.  At the moment,
C<_Householder_simple> works and is the result against which C<_Householder_fast>
should be measured.  At the moment, the only major difference comes in
how C<alpha> and C<$v> are calculated - in this version the code is more
compact.
		},
		Code => q{
			double alpha, beta, gamma;
			int j, k;
			int M = $SIZE(m);
			int N = $SIZE(n);
			threadloop %{
			loop(i) %{
				/* Note n is the number of columns, so n < m for this to work */
				loop(n) %{
					/* Compute alpha and build v */
					alpha = 0;
					for (j = n; j < M; j++) {
						alpha += $A(m => j) * $A(m => j);
						$v(m => j) = $A(m => j);
					}
					alpha = sqrt(alpha);
					if ($A(m => n) > 0)
						alpha = -alpha;
					$v(m => n) -= alpha;
					
					beta = 0;
					for (j = n; j < M; j++) {
						beta += $v(m => j) * $v(m => j);
					}
					
					if (beta != 0) {
						// apply the Householder transformation to the
						// remaining submatrix, column j, row k
						for (j = n; j < N; j++) {
							// Get the dot product of v and the jth
							// column of the submatrix, called gamma:
							gamma = 0;
							for (k = n; k < M; k++)
								gamma += $v(m => k) * $A(n => j, m => k);
							// now subtract a scaled version of v, based
							// on gamma, from the submatrix's jth
							// column:
							for (k = n; k < M; k++)
								$A(n => j, m => k) -= 2 * gamma / beta * $v(m => k);
						}
						
						// apply the Householder transformation to y
						// Get the dot product of v and y:
						gamma = 0;
						for (k = n; k < M; k++)
							gamma += $v(m => k) * $y(m => k);
						
						// subtract a scaled copy of v from y
						for (k = n; k < M; k++) {
							$y(m => k) -= (2 * gamma / beta) * $v(m => k);
						}
					}
				%}
			%}
			%}
		},
	);
};

# Note that this modifies A and y!
use Inline
	Pdlpp => q{
	pp_def('_Householder_simple',
		Pars => 'y(i, m); A(n, m); [t]v(m);',
		Code => q{
			double alpha, beta, gamma;
			int j, k;
			int M = $SIZE(m);
			int N = $SIZE(n);
			threadloop %{
			loop(i) %{
				/* Note n is the number of columns, so n < m for this to work */
				/* Loop through each column */
				loop(n) %{
					/* Compute alpha */
					alpha = 0;
					for (j = n; j < M; j++)
						alpha += $A(m => j) * $A(m => j);

					alpha = sqrt(alpha);
					if ($A(m => n) > 0)
						alpha = -alpha;

					/* Build v */
					loop (m) %{
						$v() = (m < n) ? 0 : $A();
					%}
					$v(m => n) -= alpha;
					
					/* Compute beta */
					beta = 0;
					loop (m) %{
						beta += $v() * $v();
					%}
					
					if (beta != 0) {
						// apply the Householder transformation to the
						// remaining submatrix, column j, row k
						for (j = n; j < N; j++) {
							// Get the dot product of v and the jth
							// column of the submatrix, called gamma:
							gamma = 0;
							for (k = n; k < M; k++)
								gamma += $v(m => k) * $A(n => j, m => k);
							// now subtract a scaled version of v, based
							// on gamma, from the submatrix's jth
							// column:
							for (k = n; k < M; k++)
								$A(n => j, m => k) -= 2 * gamma / beta * $v(m => k);
						}
						
						// apply the Householder transformation to y
						// Get the dot product of v and y:
						gamma = 0;
						for (k = n; k < M; k++)
							gamma += $v(m => k) * $y(m => k);
						
						// subtract a scaled copy of v from y
						for (k = n; k < M; k++) {
							$y(m => k) -= (2 * gamma / beta) * $v(m => k);
						}
					}
				%}
			%}
			%}
		},
	);
};

use Inline
	Pdlpp => q{
	pp_def('_backsub',
		Pars => 'y(i, m); A(n, m); [o]coefs(i, n);',
		Code => q{
			int j, k;
			double tmp;
			int N = $SIZE(n);
			loop (i) %{
				// Remember, n is the number of columns and must be
				// less than or equal to m.  Start at the lowest-right
				// diagonal entry:
				for (j = N-1; j > -1; j--) {
					tmp = $y(m => j);
					for (k = N-1; k > j; k--) {
						tmp -= $A(n => k, m => j) * $coefs(n => k);
					}
					$coefs(n => j) = tmp / $A(n => j, m => j);
				}
			%}
		},
	);
};


sub SVD {
	
}

# The 'curry' functions are curried (see Higher Order Perl) functions
# that I know I can trust because the logic makes sense to me.  The
# 'build' functions also return functions, but create them by evaling
# a hand-crafted string.  Presumably, the 'build' functions are faster
# than the 'curry' functions, though I have not (yet) run benchmarks
# to be sure that's true.
sub curry_anon_poly {
	my @coefs = @_;
	if (@coefs == 1 and ref ($coefs[0]) =~ /PDL/) {
		@coefs = $coefs[0]->list;
	}
	
	return sub {
		my $t = topdl(@_);
		my $result = 0;
		for (my $i = 0; $i < @coefs; $i++) {
			$result += $t**$i * $coefs[$i];
		}
		return $result;
	};
}

sub curry_anon_func {
	my ($coefs, $subs) = @_;
	my @coefs = @$coefs;
	my @subs = @$subs;
	
	# Make sure the number of subs and coefficients match
	@subs == @coefs
		or croak("The number of coefficients and functions must match!");
	
	return sub {
		my $t = topdl(@_);
		my $result = 0;
		for (my $i = 0; $i < @subs; $i++) {
			my $func = $subs[$i];
			$result += $coefs[$i] * &$func($t);
		}
		return $result;
	};
}

sub build_anon_poly {
	my @coefs = @_;
	if (@coefs == 1 and ref ($coefs[0]) =~ /PDL/) {
		@coefs = $coefs[0]->list;
	}
	my $degree = scalar(@coefs) - 1;

	# Build the function-delcaration part of the string
	my $sub_string = 'sub {
		my $t = topdl(@_);
		return (';
	
	# Build the sum, starting with the constant term:
	$sub_string .= shift @coefs;
	# If we have more terms, include them:
	while (my $coef = shift @coefs) {
		$sub_string .= "  +  \$t * ($coef";
	}
	# Add trailing parentheses and close the return statement
	$sub_string .= ')' x ($degree + 1) . ';
	}';
	
	print "Defined function:\n$sub_string\n" if $DEBUG;
	
	return eval $sub_string;
}

sub build_anon_func {
	my ($coefs, $subs) = @_;
	my @coefs = @$coefs;
	my @subs = @$subs;
	
	# Make sure the number of subs and coefficients match
	@subs == @coefs
		or croak("The number of coefficients and functions must match!");
	
	# Build the function-delcaration part of the string
	my $sub_string = 'sub {
		my $t = topdl(@_);
		return (';
	
	$sub_string .= "$coefs[0] * \&{\$subs[0]}(\$t)";
	for (my $i = 1; $i < @subs; $i++) {
		$sub_string .= "  +  $coefs[$i] * \&{\$subs[$i]}(\$t)";
	}

	# Close the anonymous sub
	$sub_string .= ")\n}";
	
	print "Defined function \n$sub_string\n" if $DEBUG;

	return eval $sub_string;
}


sub solve {
	my $arg1 = shift;
	my $arg2 = shift;


	# Here's some argument checking code, used for anonymous function building.

	my @coefs; my @subs;
	# Unpack @_, i.e. process the arguments.
	foreach (@_) {
		if (ref =~ /PDL/) {
			# If it's a piddle, assume it's a collection of coefficients
			push @coefs, $_->list;
		}
		elsif (ref =~ /CODE/) {
			# If it's a code ref, assume it's one of the basis functions
			push @subs, $_;
		}
		elsif (not ref) {
			# If it's a scalar, assume it's a coefficient
			push @coefs, $_;
		}
		else {
			croak("You should only send coefficients and subroutine references to curry_anon_sub");
		}
	}


}




1; # Magic true value required at end of module
__END__












=pod

 # In scalar context, you get the coefficients:
 my $coefs = solve(...);
 # In list context, you can get a hash with lots of goodies:
 my %results = solve(...);
 my $coefs = $results{coefs};
 my $residual = $results{residual};

 # The hash will have other items depending on the method you used.
 my $method = $results{method};  # This will either be Givens or Householder
 my $func = $results{func};  # The resulting linear combination of functions or
                             # polynomials
 
 # When solving for a linear combination of functions, you need to supply a
 # dependent variable (here called $x) that is the same size as $y.
 # You can supply a collection of functions, so this will try to find a linear
 # combination of f1() and f2() that seem to describe $y:
 my $coefs = solve($y, $x, \&f1, \&f2)
 # You can indicate a polynomial degree, so this will try to find a linear
 # combination of 1, x, and x**2 that seem to describe $y:
 my $coefs = solve($y, $x, 2);
 # For both of these forms, one of the results returned in list context is a
 # function corresponding to the fit results:
 my $fit_func = $results{function};
 
 # Sometimes you simply have a matrix that you have constructed for which you
 # want to obtain the linear least squares coefficients.  In that case, you
 # simply need to be sure that $A is one dimension bigger than $y, and the
 # number of coefficients will be determined by the size of the last dimension
 # of $A:
 my $coefs = solve($y, $A);
 
 # Finally, you can pass options into the solver if you want to fine-tune it,
 # by passing them in an anonymous hash:
 my $coefs = solve($y, $A, {method => 'Householder'});
 # At the moment, the only option is method, which takes values of Householder,
 # Givens (not yet implemented), SVD (not yet implemented).
 
 # We want to find the linear combination of functions f1, f2, etc that describe
 # our data, which I will address as x(t), so
    x(t) =~ c1 * f1(t) + c2 * f2(t) + ...
 If we have two functions and five sample times, this would be written as
 
 [ f1(t1) f2(t1) ]           [ x(t1) ]
 [ f1(t2) f2(t2) ]   [c1]    [ x(t2) ]
 [ f1(t3) f2(t3) ] * [  ] =~ [ x(t3) ]
 [ f1(t4) f2(t4) ]   [c2]    [ x(t4) ]
 [ f1(t5) f2(t5) ]           [ x(t5) ]

=cut




















=head1 NAME

PDL::Fit::Solver - A pure PDL module for performing linear data fitting using
Householder and Givens transformations, and SVD.

=head1 VERSION

This document describes PDL::Fit::Solver version 0.0.1


=head1 SYNOPSIS

 use PDL::Fit::Solver;
 
 # Get the coefficients when fitting the data to a second-order polynomial:
 $coefs = solve($y, $x, 2);
 $coefs = solve($y, $x, {order => 2});
 # Force the Givens method:
 $coefs = solve($y, $x, 2, {method => 'Givens'});
 $coefs = solve($y, $x, {order => 2, method => 'Givens'});
 # Get a whole bunch of goodies:
 %results = solve($y, $x, 2);
 my ($coefs, $residual, $func, $method)
    = @results{ qw/coefs residual func method/ };
 
 # Solve for a bunch of homemade functions:
 my $coefs = solve($y, $x, \&f1, \&f2, ...);
 
 # Solve for a hand-constructed matrix:
 my $coefs = solve($y, $A);


=head1 DESCRIPTION

This provides a module for performing linear least squares fitting.  It is
written in pure PDL.  There is no PDL::PP, XS, C, or FORTRAN code involved,
except possibly as defined in other modules.

=head2 Return Context

The return value of C<solve> depends on its context.  In scalar context, it will
return a piddle containing the coefficients of the fit.  In list context, it
will return a hash with lots of useful information.  The contents of the hash
include:

=over

=item * coefs

A piddle containing the coefficients of the fit.

=item * residual

The residual of the fit, one estimate of how well the model explains the data.

=item * method

The name of the method used for the calculation.  Should be one of 'Givens',
'Householder', or 'SVD'.

=item * func

When you supply C<solve> with a collection of functions or indicate a polynomial
order, this will contain a reference to a construction of the linear combination
of the functions.

=back

=head1 FUNCTIONS

=for author to fill in:
    Write a separate section listing the public components of the modules
    interface. These normally consist of either subroutines that may be
    exported, or methods that may be called on objects belonging to the
    classes provided by the module.



In addition to C<solve>, this module defines a number of worker functions that
the C<solve> method uses and which really do all of the work.  These include:

=head2 PDL::Fit::Solver::Householder, PDL::Fit::Solver::Givens, PDL::Fit::Solver::SVD

The actual matrix-solving functions.  You can invoke these directly, although
you could just as well specify using them in the options hash.  So, for example,
these are equivalent:

 my %results = solve($y, $A, {method => 'Givens'});
 my %results = PDL::Fit::Solver::Givens($y, $A);

However, these are matrix-only methods.  They will not accept collections of
functions.

=head2 PDL::Fit::Solver::build_anon_func

Creates and returns an anonymous function as the linear combination of the 
supplied functions using the given coefficients.  The number of coefficients
and supplied functions need not match, though a warning will be issued if
$PDL::Fit::Solver::DEBUG is set to a true value.

 my $func = build_anon_func($coefs, \&f1, \&f2);

=head2 PDL::Fit::Solver::build_anon_poly


=head1 DIAGNOSTICS

=for author to fill in:
    List every single error and warning message that the module can
    generate (even the ones that will "never happen"), with a full
    explanation of each problem, one or more likely causes, and any
    suggested remedies.

=over

=item C<< Error message here, perhaps with %s placeholders >>

[Description of error here]

=item C<< Another error message here >>

[Description of error here]

[Et cetera, et cetera]

=back


=head1 DEPENDENCIES

=for author to fill in:
    A list of all the other modules that this module relies upon,
    including any restrictions on versions, and an indication whether
    the module is part of the standard Perl distribution, part of the
    module's distribution, or must be installed separately. ]

None.


=head1 INCOMPATIBILITIES

=for author to fill in:
    A list of any modules that this module cannot be used in conjunction
    with. This may be due to name conflicts in the interface, or
    competition for system or program resources, or due to internal
    limitations of Perl (for example, many modules that use source code
    filters are mutually incompatible).

None reported.


=head1 BUGS AND LIMITATIONS

=for author to fill in:
    A list of known problems with the module, together with some
    indication Whether they are likely to be fixed in an upcoming
    release. Also a list of restrictions on the features the module
    does provide: data types that cannot be handled, performance issues
    and the circumstances in which they may arise, practical
    limitations on the size of data sets, special cases that are not
    (yet) handled, etc.

No bugs have been reported.

Please report any bugs or feature requests to
C<bug-pdl-fit-solver@rt.cpan.org>, or through the web interface at
L<http://rt.cpan.org>.


=head1 AUTHOR

David Mertens  C<< <mertens2@illinois.edu> >>


=head1 LICENCE AND COPYRIGHT

Copyright (c) 2010, David Mertens C<< <mertens2@illinois.edu> >>. All rights reserved.

This module is free software; you can redistribute it and/or
modify it under the same terms as Perl itself. See L<perlartistic>.


=head1 DISCLAIMER OF WARRANTY

BECAUSE THIS SOFTWARE IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
FOR THE SOFTWARE, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN
OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
PROVIDE THE SOFTWARE "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER
EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE
ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE SOFTWARE IS WITH
YOU. SHOULD THE SOFTWARE PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL
NECESSARY SERVICING, REPAIR, OR CORRECTION.

IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
REDISTRIBUTE THE SOFTWARE AS PERMITTED BY THE ABOVE LICENCE, BE
LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL,
OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE
THE SOFTWARE (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING
RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A
FAILURE OF THE SOFTWARE TO OPERATE WITH ANY OTHER SOFTWARE), EVEN IF
SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGES.
