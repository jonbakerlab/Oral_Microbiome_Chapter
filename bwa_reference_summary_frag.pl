#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Getopt::Long qw(GetOptionsFromArray);

my @files;
my $rate_flag;
my $verbose;
my $help;
my $status = GetOptionsFromArray(\@ARGV, 
				 "files=s{,}" => \@files,
				 #"sam" => \$sam_flag,
				 "rate" => \$rate_flag,
				 "verbose" => \$verbose,
				 "help" => \$help);

my $usage = 
    "This program takes BWA mapping files (SAM/BAM) and counts the number of mapped reads for each reference sequence.\n\n" .
    "Usage:$0\n" .
    "\t-f,--file\trequired\tstring(s)\tmapping-file1,...\n" .
    "\t-r,--rate\toptional\tflag\tratio if mulitple best hit?\t(default:off)\n" .
    "\t-v,--verbose\toptional\tflag\tverbose\t(default:off)\n" .
    "\t-h,--help\toptional\tflag\tthis message\n" .
    " e.g.:$0 -f test.bam -r\n";

if ( $help || @files == 0 ) {
    print $usage;
    exit;
}


my %Counts;
my %Sbjcts;		    ## subject matches for each query sequence
my $bwa_mem = 0;
my $first_read = 1;
foreach my $file ( @files ) {
    print STDERR "Loading $file\n";
    my $fh;
    ## Ignore unmapped reads (0x4) and supplementray alignments (0x800)
    my $exe = "samtools view -h -F 0x4 $file | samtools view -h -F 0x800 - |";
    open($fh, $exe) || die("Open error:$file");

    my $line = 0;
    my ( $prev, $curr );
    while ( <$fh> ) {
	## Get reference ids
	if ( /^\@SQ\s+SN:([^\s]+)/ ) {
	    $Counts{$1} = 0; ## Initialize all references match count to be 0
	} elsif ( /\@PG/ ) {
	    $bwa_mem = 1 if /mem/;
	}
	next if /^\@/;

	$line++;		
	parseLine( $_ );
    }
   
    print STDERR "#line:$line\n" ;
}

print STDERR "#read:" . scalar keys %Sbjcts, "\n";

## For each query sequence, increment match count to reference sequences or 
## evenly distribute to reference sequences (1/N).
foreach my $q ( keys %Sbjcts  ) {
    my @sbj = @{$Sbjcts{$q}};
    foreach my $s ( @sbj ) {
	if ( ! defined $rate_flag ) {
	    $Counts{$s}++;
	} else {
	    $Counts{$s} += 1/scalar(@sbj);
	}
    }
	
} 

foreach my $r ( sort keys %Counts ) {
    print join("\t", $r, $Counts{$r}), "\n";
}

sub getScore
{
    my ( $column ) = @_;
    $column =~ /NM:i:(\d+)/;
    return $1;
}

sub getOtherSbjcts
{
    my ( $edit, $Col) = @_;
    my @Sbj;
    ##AFNN01000024,-14652,101M,0
    if ( $Col->[$#$Col] =~ /XA:Z:(.+);$/ ) {
	my @others = split /;/, $1;
	foreach my $o ( @others ) {
	    my @f = split /,/, $o;
	    next if $f[$#f] > $edit;
	    push @Sbj, $f[0];
	}
    }
    return @Sbj;
}

sub parseLine
{
    my ($line) = @_;
    my @c = split /\s+/, $line;
	
    my $qid = $c[0];
    my $ref = $c[2];
    my $ecol = $bwa_mem ? 11 : 12;
    if ( $first_read ) {
	checkFormat($c[$ecol]);
	$first_read = 0;
    }
    my $edit = getScore($c[$ecol]);

    push @{$Sbjcts{$qid}}, $ref;
    push @{$Sbjcts{$qid}}, getOtherSbjcts( $edit, \@c );
}

sub checkFormat
{
    my ( $column ) = @_;
    $column =~ /NM:i:(\d+)/;
    die "Edit distance column error ($column)" unless $column =~ /^NM/;
}
