# Uses the Householder method, but different argument passing
# scheme, and allows for weights.

# Use exporter stuff eventually.
# Working here - document this and write tests

use strict;
use warnings;
use Carp qw(carp);
use PDL::Core ':Internal';

package PDL::Fit::OO;

# Builds a Fit object
sub new {
	my $class = shift;
	return bless shift;
}

sub get_coefs {
	return $_[0]->{coefs};
}

sub get_fit_func {
	my $this = shift;
	return sub { $this->eval_at(@_); }
}

sub get_method {
	my $this = shift;
	return $this->method;
}

sub eval_at {
	# Arguments - the object and the desired times
	my $this = shift;
	my $time = shift;

	# piddlify the incoming arg and set up the answer piddle
	my $m_time = topdl($time);
	my $to_return = zeroes($m_time);

	# Extract some useful values from the fit object
	my $coefs = $this->coefs;
	my $t = $this->{timestamp};
	
	my $i = 0;
	foreach (@{$this->{system}}) {
		if (ref eq '') {
			# It's a polynomial power, so evaluate 
			# current power (and increment that power)
			$to_return += $coefs->at($i) * $m_time**$_;
		}
		elsif (ref =~ /PDL/) {
			# It's a piddle, so interpolate
			my ($to_add) = $m_time->interpolate($t, $_);
			$to_return += $to_add * $coefs->at($i);
		}
		elsif (ref =~ /CODE/) {
			# It's a function, so evaluate it:
			$to_return += &$_($m_time) * $coefs->at($i);
		}
		$i++;
	}
	
	# Return a real scalar if they gave us a scalar
	return $to_return->at(0) if (ref($time) eq '');
	# Otherwise return the piddle
	return $to_return;
}



use PDL::Fit::Householder;

sub PDL::fit {
	my $data = shift;
	
	# Get the args
	my ($stuff_to_fit, $t, $weights) = _process_args(@_);
	
	# Make sure they passed something to fit
	barf("You didn't pass anything to fit!") unless @$stuff_to_fit;

	# Make sure that if functions or polynomial fits are being performed
	# that we have a dependent variable
	if (grep {ref eq '' or ref eq 'CODE'} @$stuff_to_fit) {
		barf("You must specify a dependent variable for function and polynomial fitting")
			unless defined $t;
	}
	
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
	my $A = cat(@piddles_to_fit)->transpose * $weights->transpose;
	my $y = $data->copy->dummy(0) * $weights->transpose;

	# Perform the fit
	$y->_Householder($A);
	my $coefs = $y->_backsub($A)->squeeze;
	
	# Eventually return coefficients, residual, method, etc if requested
	my ($residual);
	return PDL::Fit::OO->new(
			coefs => $coefs,
			timestamp => $t,
			system => $stuff_to_fit,
			residual => $residual,
			method => 'Householder',
			diagonalized_matrix => $A,
			manipulated_results => $y,
		);
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
				eval { my $a = sequence(2,3,4) + &$func(sequence(2,3,4)) };
				barf("Function $i generates output that is not compatible with input")
					if ($@);
			}
			
			# If we're here then the list is good
			push @stuff_to_fit, @funcs;
			
			next ARGUMENT;
		}
		
		# Check for fit-system, in which case we check that it's an array
		# ref but otherwise assume it's correct and let it override
		# anything else.
		elsif (/^fit_system$/) {
			barf("fit systems must be array references")
				unless ref ($value) eq 'ARRAY';

			@stuff_to_fit = @$value;
			last ARGUMENT;
		}
		
		# Warn on unknown options
		else {
			carp ("PDL::Fit: Unknown option $_; ignoring");
		}
	}

	return (\@stuff_to_fit, $t, $weights);
}

1;

=pod

This needs documentation. In the meantime, this will have to suffice.

=head1 SYNOPSIS

Here's an example that uses a standard fit:

 use PDL;
 use PDL::Fit::OO;
 
 my $t = sequence(100)/10;
 my $signal = 4 * sin($t) - 3 * cos($t) + 0.1 * grandom($t) + 2;
 
 # fit the data to a cosine, sine, and linear polynomial
 my $fit = $signal->fit(functions => [\&PDL::sin, \&PDL::cos]
                         , poly => 1, t => $t);
 
 print $fit->get_coefs;
 # Should look something like 4, -3, 2, 0
 
 # Plot the results for comparison
 use aliased 'PDL::Graphics::PLplot';
 my $pl = PLplot->new(DEV => YOUR_FAVORITE_DEVICE_HERE)
 # Plot the data
 $pl->xyplot($t, $signal, PLOTTYPE => 'POINTS');
 # Plot the fit
 $pl->xyplot($t, $fit->eval_at($t), PLOTTYPE => 'LINE', COLOR => 'RED');
 

Here's an example that uses weighted fitting:

 use strict;
 use warnings;
 use PDL;
 use PDL::Fit::OO;
 
 my $t = sequence(60)/10;
 my $signal = cos($t)->remember_fit_system;
 
 # First try the whole thing, evenly weighted
 my $unweighted_fit = $signal->fit(polynomial => 2, t => $t);
 
 # Next try to focus on the center, new pi
 my $weights = exp(-($t - 3.14)**2);
 my $weighted_fit = $signal->fit(
                 polynomial => 2
               , t => $t
               , weights => $weights
           );
 
 # Plot the results
 use aliased 'PDL::Graphics::PLplot';
 my $pl = PLplot->new(DEV => 'xwin');
 # Plot the data
 $pl->xyplot($t, $signal, PLOTTYPE => 'POINTS');
 print "Plotting the unweighted fit in red\n";
 # Plot the unweighted fit in red
 $pl->xyplot($t, $unweighted_fit->eval_at($t), PLOTTYPE => 'LINE', COLOR => 'RED');
 print "Plotting the weighted fit in blue\n";
 # Plot the weighted fit in blue
 $pl->xyplot($t, $weighted_fit->eval_at($t), PLOTTYPE => 'LINE', COLOR => 'BLUE');

=cut
