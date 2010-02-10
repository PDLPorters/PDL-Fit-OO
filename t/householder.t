#!perl

use Test::More tests => 1;
use PDL::Fit::Householder;

# The following tests need to happen:
# 1) check argument parsing for sub Householder
# 2) check _check_piddle to make sure it works right
# 3) check that _backsub works correctly by constructing a matrix and y,
#    finding coefs, and then making sure that $matrix x $coefs == $y
# 4) check that _householder gives the same results as PDL::MatrixOps
#    for diagonal matrices
# 5) Check that generated signals with noise can be fit with reasonable
#    results.
