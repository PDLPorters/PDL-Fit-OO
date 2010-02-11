use strict;
use warnings;

use Module::Build::PDL;

my $builder = Module::Build::PDL->new
(
	config_requires => {
		'Module::Build::PDL'  => 0,
	},
    module_name         => 'PDL::Fit::Householder',
    license             => 'perl',
    dist_author         => 'David Mertens <mertens2@illinois.edu>',
    dist_version_from   => 'lib/PDL/Fit/Householder.pd',
    requires => {
        'version'    => 0,
        'PDL'        => 0,
    },
    build_requires => {
        'Test::More'       => 0,
    	'Test::Exception'  => 0,
    	'Test::Warn'       => 0,
    	'Test::LectroTest' => 0,
    },
    add_to_cleanup      => [ 'PDL-Fit-Householder-*' ],

);
$builder->create_build_script;
