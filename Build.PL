use strict;
use warnings;

#use lib '.';
use Module::Build::PDL;

my $builder = Module::Build::PDL->new
(
    module_name         => 'PDL::Fit::Solver',
    license             => 'perl',
    dist_author         => 'David Mertens <mertens2@illinois.edu>',
    dist_version_from   => 'lib/PDL/Fit/Solver.pm',
    requires => {
        'Test::More' => 0,
        'version'    => 0,
    },
    add_to_cleanup      => [ 'PDL-Fit-Solver-*' ],

);
$builder->create_build_script;
