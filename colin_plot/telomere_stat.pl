#!/usr/bin/perl

use strict;
use warnings;

my $genome_file = $ARGV[0];
my $genome = load_seq($genome_file);

my $telo_string = 'CCCTAACCCTAACCCTAA';
my $revcom_telo_string = revcom_seq($telo_string);
my $telo_repeat_len = length $telo_string;
my $max_telo_dist = 50000;

my %telo_repeat = ();
foreach my $chr (sort by_strnum_td keys $genome) {
	my $chr_len = length $genome->{$chr};
	my $max_len = $chr_len - $telo_repeat_len - 1;
	foreach my $i (0 .. $max_len) {
		my $str = substr($genome->{$chr}, $i, $telo_repeat_len);
		if ($str eq $telo_string or $str eq $revcom_telo_string) {
			if ($i < $max_len/2 and $i < $max_telo_dist) {
				$telo_repeat{$chr}{'left'} = $i;
			}
			if ($i > $max_len/2 and $max_len - $i < $max_telo_dist) {
				$telo_repeat{$chr}{'right'} = $i;
			}
		}
	}
	$telo_repeat{$chr}{'chrlen'} = $chr_len;
	if (!exists $telo_repeat{$chr}{'left'}) {
		$telo_repeat{$chr}{'left'} = 'nd';
	}
	if (!exists $telo_repeat{$chr}{'right'}) {
		$telo_repeat{$chr}{'right'} = 'nd';
	}
}

print "ref\tchr_len\tleft\tright\n";
foreach my $chr (sort by_strnum_td keys $genome) {
	my $num = $chr;
	$num =~ s/\S+\_//;
	print $num, "\t", $telo_repeat{$chr}{'chrlen'}, "\t", $telo_repeat{$chr}{'left'}, "\t", $telo_repeat{$chr}{'right'}, "\n";
}


sub load_seq {
	my $seq = shift;
	my %seq = ();
	open (IN, $seq) or die "Cannot open file $seq: $!\n";
	my $seq_id = undef;
	while (<IN>) {
		chomp;
		if (/^\>(\S+)/) {
			$seq_id = $1;
		} else {
			s/\s//g;
			$seq{$seq_id} .= uc($_);
		}
	}
	close IN;
	return \%seq;
}

sub revcom_seq {
	my ($seq) = @_;
	$seq = reverse $seq;
	$seq =~ tr/ATCGNXatcgnx/TAGCNXtagcnx/;
	return $seq;
}

sub by_strnum_td {
	$a =~ /\_(\d+)/;
	my $numa = $1;
	$b =~ /\_(\d+)/;
	my $numb = $1;
	return $numa <=> $numb;
}

sub by_num {
	$a <=> $b;
}
