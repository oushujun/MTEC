#!/usr/bin/perl -w
use strict;
#Shujun Ou (shujun.ou.1@gmail.com; 11/22/2019)

my $usage = "
	Rename short TIR elements to MITE. Only TIRs are renamed if other class exists\n
	perl rename_MITE.pl TE.fasta > TE.fasta.renamed\n\n";

my $len_cutoff = 600; #TIR elements <= this value are classified as MITEs, others DNA TEs.

$/ = "\n>";
while (<>){
	s/>//g;
	my ($id, $seq) = (split /\n/, $_, 2);
	$seq =~ s/\s+//g;
	my ($name, $class, $superfam) = ($1, $2, $3) if $id =~ /^(.*)#(.*)\/(.*)$/;
	#ZM00001_consensus#DNA/DTA
	if (length $seq < $len_cutoff and $class =~ /DNA|MITE/ and $superfam !~ /Helitron|DHH/){
		$class = "MITE";
		}
#print "$id\t$name, $class, $superfam\n";
	print ">$name#$class/$superfam\n$seq\n";
	}

