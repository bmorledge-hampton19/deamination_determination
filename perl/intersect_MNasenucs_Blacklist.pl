#!/usr/bin/perl 
# intersect_TFBS_Blacklist_100bp_window.pl

use strict;
use warnings;

print STDERR "Enter filename of TFBS to intersect with BLIST coordinates\n";
my $tfbs_file = <STDIN>;
chomp $tfbs_file;
open ( TFBS, $tfbs_file ) || die "Couldn't open file: $tfbs_file\n";

# open BLIST coordinates
print STDERR "Enter filename of BlackList coordinate file:\n";
my $blist_file = <STDIN>;
chomp $blist_file;
open (BLIST, $blist_file) || die "Couldn't open promoter file: $blist_file\n";
print STDERR "processing BlackList coordinates...\n";

my $window = 73;
# process promoters into lookup hash (1 => promoter, undef/not exist ==> not promoter)
my %BLIST_coord;
my $chr = "";
while ( my $line = <BLIST> )
{
	chomp $line;
	my @fields = split /\t/, $line;
	if ( $fields[0] =~ /(chr[0-9A-Za-z_]+)/ )
	{
		my $temp = $1;
		if ( $chr ne $temp )
		{
			print STDERR "Starting to process $temp\n";
			$chr = $temp;
		}

		push @{$BLIST_coord{$chr}}, [$fields[1], $fields[2]];
	}
	else
	{
		print STDERR "Misformatted line: $line\n";
	}
}
close (BLIST);

# intersect TF coordinates with BLIST --> any overlap means active BLIST TF

$tfbs_file =~ s/\.txt//;

my $InBLIST_file = $tfbs_file . "_BLIST.txt";
my $notInBLIST_file = $tfbs_file . "_notBLIST.txt";

open (BLIST, ">$InBLIST_file") || die "Couldn't open InBLIST file: $InBLIST_file\n";
open (NOTBLIST, ">$notInBLIST_file") || die "Couldn't open NOTInBLIST file: $notInBLIST_file\n";

print STDERR "Starting processing TFBS coordinates...\n";
$chr = "";
my %BLIST_lookup;
while ( my $line = <TFBS>)
{
	chomp $line;
	my @fields = split /\t/, $line;
        if ( $fields[0] =~ /(chr[0-9A-Za-z_]+)/ )
        {
		my $temp = $1;
		if ( $temp ne $chr )
		{
                	print STDERR "Starting to process $temp\n";
			$chr = $temp;
			# empty hash
			%BLIST_lookup = ();
			foreach my $blist ( @{$BLIST_coord{$chr}} )
			{
				my $start = $blist->[0];
				# start coordinate is 0-based, so convert to 1-based for consistency
				$start++;
				my $end = $blist->[1];
                		for ( my $i = $start; $i <= $end; $i++ )
                		{
                        		# Define this nucleotide as blacklisted 
                        		$BLIST_lookup{$fields[0]}{$i} = 1;
                		}
			}
		}
                my $nucdyad = $fields[2]; # already 1-based 
		if ( $nucdyad - 1 != $fields[1] )
		{	die "Weird dyad: $temp\n";	}

		# expand size of window
		my $winstart = $nucdyad - $window;
		my $winend = $nucdyad + $window;		

		my $prox_flag = 0;
		for ( my $i = $winstart; $i <= $winend; $i++ )
		{
			if ( exists $BLIST_lookup{$chr}{$i} && $BLIST_lookup{$chr}{$i} == 1 )
			{
				$prox_flag = 1;	
			}
		}

                if ( $prox_flag )
                {
                        print BLIST "$line\n";
                }
                else
                {
                        print NOTBLIST "$line\n";
                }
	}
	else 
	{
		die "Misformatted line: $line\n";
	}
}
