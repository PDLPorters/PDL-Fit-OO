package PDL::Demos::Householder;
# A bunch of fitting demos
# Copyright (C) 2010 David Mertens
# Permission is granted to distribute this demo as long as
# as proper atribution to the author is given.

use PDL;
use PDL::Fit::Householder;

# Need to run these in combination since PDL::Demos::Screen holds
# the package PDL::Demos::Routines!
use PDL::Demos::Screen;
PDL::Demos::Routines->import();

sub comment($);
sub act($);
sub output;

sub run {

comment <<'INTRO';
  This is a tour of the fitting capabilities of PDL::Fit::Householder.
  Hopefully it will give you a good idea of how to use this module to
  fit your own data.

  As with all modules, you should start by including it:
      use PDL;
      use PDL::Fit::Householder;

  In the examples that follow, I will often resort to manufacturing my
  data.  In the real world, of course, you would use Householder to fit
  data you have from other sources.

INTRO

act <<'POLYNOMIAL';
  # We'll begin with a polynomial fit.  First a dependent variable:
  $t = sequence(1000) / 50;
  # Next some mocked data:
  $data = -0.5 + 2 * $t + 5 * $t**2 + grandom($t);
  
  # This is the command to fit a second-order polynomial to $data using
  # dependent variable $t:
  $coefs = $data->Householder(2, $t);
  print $coefs, "\n";

POLYNOMIAL

act <<'FUNCTION';
  # We'll begin with a polynomial fit.  First a dependent variable:
  $t = sequence(1000) / 50;
  # Next some mocked data:
  $data = -0.5 + 2 * $t + 5 * $t**2 + grandom($t);
  
  # This is the command to fit a second-order polynomial to $data using
  # dependent variable $t:
  $coefs = $data->Householder(2, $t);
  print $coefs, "\n";

FUNCTION

comment <<'FINISH_WITH_T';
  

FINISH_WITH_T

comment <<'CONCLUSION';
  Thanks for following allong, and I hope that gives you some
  idea of how to use PDL::Fit::Householder.
CONCLUSION

}

# If we're not running under the Perldl shell, then run the demo:
run if (not defined $PERLDL::PROMPT);

1;
