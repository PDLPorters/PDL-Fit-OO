use strict;
use warnings;

use Module::Build::PDL;

my $builder = Module::Build::PDL->new
(
    module_name         => 'PDL::Fit::Householder',
    license             => 'perl',
    dist_author         => 'David Mertens <mertens2@illinois.edu>',
    dist_version_from   => 'lib/PDL/Fit/Householder.pd',
    requires => {
        'Test::More' => 0,
        'version'    => 0,
        'PDL'        => 0,
    },
    add_to_cleanup      => [ 'PDL-Fit-Householder-*' ],

);
$builder->create_build_script;
