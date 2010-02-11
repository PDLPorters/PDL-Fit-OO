#!perl

eval {require Test::Perl::Critic};

if ($@) {
	use Test::More skip_all => "Test::Perl::Critic required for testing PBP compliance, but couldn't find it";
}

Test::Perl::Critic::all_critic_ok();
