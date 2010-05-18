# Uses the Householder method, but different argument passing
# scheme, and allows for weights.

# Use exporter stuff eventually.

use PDL;
use PDL::Fit::Householder;
use strict;
use warnings;
use Carp qw(carp);
use PDL::Core ':Internal';

sub PDL::fit {
	my $data = shift;
	
	# Get the args and be sure we have an even number of them.
	my ($stuff_to_fit, $t, $weights) = _process_args(@_);
	
	# Set default weights to 1
	$weights = pdl(1) unless defined $weights;
	
	# Build the weighted matrix and copy data for fitting
	my @piddles_to_fit;
	foreach(@$stuff_to_fit) {
		if (not ref) {
			push @piddles_to_fit, $t ** $_;
		}
		elsif (ref ($_) =~ /CODE/) {
			push @piddles_to_fit, &$_($t);
		}
		elsif (ref ($_) =~ /PDL/) {
			push @piddles_to_fit, $_;
		}
	}
	my $A = cat(@piddles_to_fit)->transpose / $weights->transpose;
	my $y = $data->copy->dummy(0) / $weights;

	# Perform the fit
	$y->_Householder($A);
	my $coefs = $y->_backsub($A)->squeeze;
	
	# Eventually return coefficients, residual, method, etc if requested
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

sub _process_args {
	my @args = @_;
	
	barf("Arguments to PDL::fit must be in key/value pairs\n")
		unless @args % 2 == 0;
	
	# These keep track of ...
	my @stuff_to_fit;		# everything we need to fit
	my $t;					# dependent variable
	my $weights;			# fitting weights
	my $dim_check_piddle;	# a piddle against which we check dimensions
							#		of other piddles

	# Begin by checking if they specified a dependent variable.
	# If so, pull it out and store it.
	TCHECK: for (my $i = 0; $i < @args; $i += 2) {
		my ($key, $value) = @args[$i, $i+1];
		
		# OK, this ain't pretty, but it's backward compatible as
		# hell and everybody knows what's going on:
		local $_ = $key;
		
		if (/^t$/ || /^dep/) {
			splice (@args, $i, 2);
			$t = $value;
			
			# working here - check that $t is a piddle; improve this
			ref ($t) =~ /PDL/
				or barf("Dependent variable is not a piddle!");
			
			$dim_check_piddle = $t;
			
			last TCHECK;
		}
	}
	
	# Next, pull out the weights, if given
	WEIGHTCHECK: for (my $i = 0; $i < @args; $i += 2) {
		my ($key, $value) = @args[$i, $i+1];
		
		# OK, this ain't pretty, but it's backward compatible as
		# hell and everybody knows what's going on:
		local $_ = $key;
		
		if (/^weight/) {
			splice (@args, $i, 2);
			$weights = $value;
			
			# working here - check that $weights is a piddle
			# improve this
			ref ($weights) =~ /PDL/
				or barf("Weights variable is not a piddle!");
			
			# Ensure compatibility with $t, if $t is defined
			if (defined $dim_check_piddle) {
				eval { my $a = $dim_check_piddle + $weights };
				barf("Weights and dependent variable have incompatible dimensions")
					if $@;
			}
			# and set the reference piddle to weights if $t is not defined
			else {
				$dim_check_piddle = $weights;
			}
			
			last WEIGHTCHECK;
		}
	}
	
	
	# Run through the key/value pairs
	ARGUMENT: for(my $i = 0; $i < @args; $i += 2) {
		my ($key, $value) = @args[$i, $i+1];
		local $_ = $key;

		# Check for polynomial order spec:
		if (/^poly/) {
			# Avoid piddle/scalar mixups
			$value = $value->at(0) if ref ($value) =~ /PDL/;
			
			# Check that the order is a non-negative integer
			barf("Polynomial order $value may not be a non-negative integer")
				unless $value == int $value and $value >= 0;
			
			# Check that we have a dependent variable specified
			barf("Fitting to a polynomial requires a dependent variable")
				unless defined $t;

			push @stuff_to_fit, $_ foreach (0 .. $value);
			
			next ARGUMENT;
		}
		
		# Check for piddles spec:
		elsif (/^piddle/) {
			my @piddles;
			if (ref ($value) =~ /PDL/) {
				@piddles = ($value);
			}
			elsif (ref ($value) eq 'ARRAY') {
				@piddles = @$value;
			}
			else {
				barf("Piddle arguments must be a piddle or an array of piddles");
			}
			
			# Make sure they actually gave us some piddles
			@piddles > 0 or barf("You didn't actually send me any piddles");
			
			# Check each piddle in the list
			for (my $i = 0; $i < @piddles; $i++) {
				my $piddle = $piddles[$i];
				
				# Make sure it's actually a piddle
				if (ref ($piddle) !~ /PDL/) {
					barf("Piddle number $i is not actually a piddle");
				}
				
				# Working here - check that piddle dimensions are
				# compatible with each other and dependent variable
				if ($i > 0) {
					eval { my $a = $piddle + $piddles[$i-1] };
					barf("Piddle $i is incompatible with neighboring piddle")
						if ($@);
				}
				elsif (defined $dim_check_piddle) {
					eval { my $a = $piddle + $dim_check_piddle };
					barf("Piddle 0 is incompatible with dependent variable or weights")
						if ($@);
				}
			}
			
			# If we're here then the list is good
			push @stuff_to_fit, @piddles;
			
			next ARGUMENT;
		}
		
		# Check for function references
		elsif (/^func/) {
			# Make sure we have a dependent variale
			barf("Fitting to functions requires a dependent variable")
				unless defined $t;
			
			my @funcs;
			if (ref ($value) eq 'CODE') {
				@funcs = ($value);
			}
			elsif (ref ($value) eq 'ARRAY') {
				@funcs = @$value;
			}
			else {
				barf("Function arguments must be an individual or array of function references");
			}
			
			# Make sure they actually gave us something
			@funcs > 0 or barf("You didn't actually send me any functions");
			
			# Check each function in the list
			for (my $i = 0; $i < @funcs; $i++) {
				my $func = $funcs[$i];
				
				# Make sure it's actually a function
				unless (ref ($func) eq 'CODE') {
					barf("Function $i is not actually a function reference");
				}
				
				# Check that the function gives results with dimensions
				# that are compatible with $t
				eval { my $a = $dim_check_piddle + &$func($t) };
				barf("Function $i outputs a piddle with dimensions that are incompatible with the dependent variale or weights")
					if ($@);
			}
			
			# If we're here then the list is good
			push @stuff_to_fit, @funcs;
			
			next ARGUMENT;
		}
		
		# Warn on unknown options
		else {
			carp ("PDL::Fi: Unknown option $_; ignoring");
		}
	}
	
	barf("You didn't pass anything to fit!") unless @stuff_to_fit;

	return (\@stuff_to_fit, $t, $weights);
}

sub PDL::function_from_fit {
	my $coefs = shift;

	# Get the args
	my ($stuff_to_fit, $t) = _process_args(@_);
	
	# If they passed in piddles for me to interpolate, I need the
	# dependent variable
	if (grep {ref =~ /PDL/} @$stuff_to_fit) {
		barf("You must specify a dependent variable for function_from_fit")
			unless defined $t;
	}
	
	# Return a function that does the interpolation:
	return sub {
		my $time = shift;
		my $poly_power = 0;
		my $to_return = zeroes(topdl($time));
		
		my $i = 0;
		foreach (@$stuff_to_fit) {
			if (ref eq '') {
				# It's a polynomial power, so evaluate 
				# current power (and increment that power)
				$to_return += $coefs->at($i) * $time**$_;
			}
			elsif (ref =~ /PDL/) {
				# It's a piddle, so interpolate
				my ($to_add) = $time->interpolate($t, $_);
				$to_return += $to_add * $coefs->at($i);
			}
			elsif (ref =~ /CODE/) {
				# It's a function, so evaluate it:
				$to_return += &$_($time) * $coefs->at($i);
			}
		}
		
		# Return a real scalar if they gave us a scalar
		return $to_return->at(0) if (ref($time) eq '');
		# Otherwise return the piddle
		return $to_return;
	};
}

=pod

  Householder(poly_order => SCALAR,
    piddle => PIDDLE | piddles => [PIDDLE, PIDDLE, ...]
    function => FUNCREF | functions => [FUNCREF, FUNCREF, ...],
    dep_var => PIDDLE | t => PIDDLE)
