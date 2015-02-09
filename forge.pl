#!/usr/bin/env perl

=head1 NAME

forge.pl - Functional element Overlap analysis of the Results of GWAS Experiments.

=head1 SYNOPSIS

forge.pl options (-f file) (-snps snplist)

=head1 DESCRIPTION

Analyse a set of SNPs for their overlap with DNase 1 hotspots compared to matched background SNPs. Identifies enrichment in DHS by tissue and plots graphs and table to display. Arbitrarily a minumum of 20 SNPs is required.  Note that if no SNPs are given the script will run on Pulmonary funciotn GWAS as an example output.

Several outputs are made.

A straight base R graphics pdf chart of the data.

A polychart (https://github.com/Polychart/polychart2) interactive javascript graphic using rCharts (http://ramnathv.github.io/rCharts/).

A dimple (http://dimplejs.org) d3 interactive graphic using rCharts.

A table using the Datatables (https://datatables.net) plug-in for the jQuery Javascript library, again accessed through rCharts.

In each of the graphics the colouring should be consistent. Blue (p > 0.05, light red or pink (0.05 >= p > 0.01), red or dark red (p <= 0.01 ) for the 95% and 99% cIs after Boonferroni correction. Or whatever other thresholds are specified.

Forge functions, plotting options and stats are provided by Forge::Forge, Forge::Plot and Forge::Stats modules.

=head1 OPTIONS

=over

=item B<data>

Data set to analyse. Either ENCODE data ('encode') or Roadmap Epigenome data ('erc'). erc by default.

=item B<peaks>

Use peaks instead of hotspots. Peaks are more stringent DNase1 peaks calls representing DNase hypersensitive sites, rather than hotspots which are regions of generalised DNAs1 sensitivity or open chromatin. Default is to use hotspots.

=item B<bkgd>

Specify whether the background matches should be picked from general set of arrays used in GWAS ('gwas') or from the Illumina_HumanOmni2.5 ('omni'). General GWAS arrays include

Affy_GeneChip_100K_Array
Affy_GeneChip_500K_Array
Affy_SNP6
HumanCNV370-Quadv3
HumanHap300v2
HumanHap550v3.0
Illumina_Cardio_Metabo
Illumina_Human1M-duoV3
Illumina_Human660W-quad

Defaults to 'gwas'. In both cases SNPs have to be on the arrays AND in the 1000 genomes phase 1 integrated call data set at phase1/analysis_results/integrated_call_sets.

=item B<label>

Supply a label that you want to use for the plotting titles, and filenames.

=item B<f>

Supply the name of a file containing a list of SNPs. Format must be given by the -format flag. If not supplied the analysis is performed either on snps provided as rsids in a comma separated list through the snps option or on a set of data from a gwas study on Pulmonary_function (http://www.ncbi.nlm.nih.gov/pubmed/21946350, http://www.ncbi.nlm.nih.gov/pubmed/20010835 and http://www.ncbi.nlm.nih.gov/pubmed/20010834). Note that 20 SNPs are required at a minimum.

=item B<snps>

Can provide the snps as rsids in a comma separated list.

=item B<min_snps>

Specify the minimum number of SNPs to be allowed. Default is 5 now are using binomial test.

=item B<thresh>

Alter the default binomial p value thresholds. Give a comma separate list of two e.g. 0.05, 0.01 for the defaults

=item B<format>

if f is specified, specify the file format as follow:

rsid = list of snps as rsids each on a separate line. Optionally can add other fields after the rsid which are ignored, unless the pvalue filter is specified, in which case Forge assumes that thesecond field is the minus log10 pvalue

bed  = File given is a bed file of locations (chr\tbeg\tend) aka Personal Genome SNP format.  bed format should be 0 based and the chromosome should be given as chrN. Hoever will also accept chomosomes as just N (ensembl) and 1-based format where beg and end are the same

vcf = File given is a vcf file.

tabix = File contains SNPs in tabix format.

ian = 1-based chr\tbeg\tend\trsid\tpval\tminuslog10pval

=item B<filter>

Set a filter on the SNPs based on the -log10 pvalue.  This works for files in the 'ian' or 'rsid' format. Give a value as the lower threshols and only SNPs with -log10 pvalues >= to the threshold will be analysed. Defaiult is no filtering.

=item B<bkgrd>

Output background stats for investigation.

=item B<reps>

THe number of background matching sets to pick and analyse. Default 1000.

=item B<ld>

Apply filter for SNPs in LD at either r2 >= 0.8 ("high LD"), or r2 >= 0.1 ("independent SNPs"). Specify ld 0.8, or ld 0.1. Default is to filter at r2 >= 0.8.  With ld filter specified, forge will report SNPs removed due to LD with another SNP in the list and will randomly pick one for each LD block.

To turn off LD filtering specify -nold

=item B<nold>

Turn off LD filtering.

=item B<noplot>

Just make the data file, don't plot.

=item B<overlap>

Just run overlaps, no enrichment analysis.

=item B<help|h|?>

Print a brief help message and exits.

=item B<man|m>

Print this perldoc and exit.

=back

=head1 LICENCE

forge.pl Functional analysis of GWAS SNPs

Copyright (C) 2013  EMBL - European Bioinformatics Institute

This program is free software: you can redistribute it and/or modify it under the terms of
the GNU General Public License as published by the Free Software Foundation, either version 3
of the License, or (at your option) any later version. This program is distributed in the hope
that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. Neither
the institution name nor the name forge.pl can be used to endorse or promote products derived from
this software without prior written permission. For written permission, please contact
dunham@ebi.ac.uk. Products derived from this software may not be called forge.pl nor may forge.pl
appear in their names without prior written permission of the developers. You should have received
a copy of the GNU General Public License along with this program.  If not, see http://www.gnu.org/licenses/.

=head1 AUTHOR

Ian Dunham, EBI

=head1 CONTACT

Ian Dunham <dunham@ebi.ac.uk>

=cut

use strict;
use 5.010;
use warnings;
use DBI;
use Sort::Naturally;
use Cwd;
use Getopt::Long;
use File::Basename;
use Config::IniFiles;
use Pod::Usage;
use Scalar::Util qw(looks_like_number);
use Forge::Stats;
use Forge::Plot;
use Forge::Forge;

my $cwd = getcwd;

my ($bkgd, $data, $peaks, $label, $file, $format, $min_snps, $bkgrdstat, $noplot, $reps, $help, $man, $thresh, $ld, $nold, $filter, $overlap, @snplist);

GetOptions (
    'data=s'     => \$data,
    'peaks'      => \$peaks,
    'bkgrd'      => \$bkgrdstat,
    'label=s'    => \$label,
    'f=s'        => \$file,
    'format=s'   => \$format,
    'snps=s'     => \@snplist,
    'min_snps=i' => \$min_snps,
    'noplot'     => \$noplot,
    'reps=i'     => \$reps,
    'thresh=s'   => \$thresh,
    'ld=f'       => \$ld,
    'nold'       => \$nold,
    'filter=f'   => \$filter,
    'overlap'    => \$overlap,
    'help|h|?'   => \$help,
    'man|m'      => \$man,

);

pod2usage(1) if ($help);
pod2usage(-verbose => 2) if ($man);

# the minimum number of snps allowed for test. Set to 5 now we have binomial p?
unless (defined $min_snps){
    $min_snps = 5;
}
if (defined $overlap) {
    $min_snps = 1;
    $nold = 1;
}

# define which data we are dealing with for the bitstrings, erc or encode
unless (defined $data ){
    $data = "erc";
}
# Label for plots
unless (defined $label){
    $label = "No label given";
}
(my $lab = $label) =~ s/\s/_/g;
$lab = "$lab.$data";
#format for reading from file
unless (defined $format){
    $format = 'rsid';
}

# Read the config file, forge.ini
my $dirname = dirname(__FILE__);
my $cfg = Config::IniFiles->new( -file => "$dirname/forge.ini" );
my $datadir = $cfg->val('Files', 'datadir');

# percentile bins for the bkgrd calculations. This is hard coded so there are enough SNPs to choose from, but could later be altered.
my $per = 10;
# number of sets to analyse for bkgrd.
unless (defined $reps){
    $reps = 1000;
}
# Which arrays to use for background
unless (defined $bkgd){
    $bkgd = "gwas";
}
# Define the thresholds to use.
my ($t1, $t2);
if (defined $thresh){
    ($t1, $t2) = split(",", $thresh);
    unless (looks_like_number($t1) && looks_like_number($t2)){
        die "You must specify numerical pvalue thresholds in a comma separated list";
    }
}
else{
    $t1 = 0.05; # set binomial p values, bonferroni is applied later based on number of samples (cells)
    $t2 = 0.01;
}

# Set r2 LD thresholds
my $r2;
unless (defined $nold){
    unless (defined $ld){
        $ld = 0.8;
    }
    unless ($ld == 0.1 || $ld == 0.8){
        die "You have specified LD filtering, but given an invalid value $ld. the format is ld 0.8, or ld 0.1";
    }
    ($r2 = $ld) =~ s /\.//;
    $r2 = "r".$r2;
}

my $dsn;
if (defined $peaks){
    $dsn = "dbi:SQLite:dbname=" . $datadir . "forge_peaks.db";
}
else{
    $dsn = "dbi:SQLite:dbname=" . $datadir . "forge.db";
}
my $dbh = DBI->connect($dsn, "", "") or die $DBI::errstr;
# snps need to come either from a file or a list
my @snps;

# A series of data file formats to accept.

if (defined $file){
    if (defined $filter) {
        unless ($format eq "ian" or $format eq "rsid"){
            warn "You have specified pvalue filtering, but this isn't implemented for files of format $format. No filtering will happen."
        }
    }
    my $sth = $dbh->prepare("SELECT rsid FROM bits WHERE location = ?");
    open my $fh, "<", $file or die "cannot open file $file : $!";
    @snps = process_file($fh, $format, $sth, $filter);
}
elsif (@snplist){
    @snps = split(/,/,join(',',@snplist));
}
else{
# Test SNPs from gwascatalog_21_03_2012  Pulmonary_function.snps.bed
# If no options are given it will run on the default set of SNPs
    warn "No SNPs given, so running for example on Pulmonary function set from the GWAS catalogue.";
    @snps = qw(rs2865531 rs2395730 rs12914385 rs11168048 rs1529672 rs357394 rs13147758 rs3769124 rs2647044 rs12504628 rs1541374 rs2869967 rs1928168 rs3094548 rs3867498 rs9978142 rs4762767 rs6903823 rs11172113 rs9310995 rs2571445 rs2070600 rs11727189 rs3734729 rs2906966 rs1036429 rs16909898 rs3995090 rs12477314 rs2544527 rs2284746 rs993925 rs2277027 rs1344555 rs1455782 rs2855812 rs2838815 rs11001819 rs12716852 rs2798641 rs4129267 rs7068966 rs12899618 rs153916 rs1551943 rs730532 rs1980057 rs3820928 rs2036527 rs10516526 rs2857595 rs3817928 rs310558 rs808225 rs12447804);
}

# Remove redundancy in the input

my %nonredundant;
foreach my $snp (@snps){
    $nonredundant{$snp}++;
}
foreach my $snp (keys %nonredundant){
    if ($nonredundant{$snp} > 1) {
        say "$snp is present " . $nonredundant{$snp} . " times in the input. Analysing only once."
    }
}
@snps = keys %nonredundant;
my @origsnps = @snps;

# Perform ld filter unless -nold is specified.

my ($ld_excluded, $output, $input);
unless (defined $nold) {
    $input = scalar @snps;
    ($ld_excluded, @snps) = ld_filter(\@snps, $r2, $dbh);
    $output = scalar @snps;
    #say join("\t", @snps);
}

# Check we have enough SNPs
if (scalar @snps < $min_snps){
    pod2usage(-verbose => 2, -message => "Fewer than $min_snps SNPs. Analysis not run\n\n", -noperldoc => 1);
}

# get the cell list array and the hash that connects the cells and tissues
my ($cells, $tissues) = get_cells($data, $dbh);

# get the bit strings for the test snps from the database file
my $rows = get_bits(\@snps, $dbh);

# unpack the bitstrings and store the overlaps by cell.
my $test = process_bits($rows, $cells, $data);

# generate stats on the background selection
if (defined $bkgrdstat){
    bkgrdstat($test, $lab, "test");
}

# Identify SNPs that weren't found and warn about them.
my @missing;
foreach my $rsid (@origsnps){
    if (defined $ld) {
        next if exists $$ld_excluded{$rsid}; # if the snps are not in the 1000 genomes set they will not be found by the LD filter, so we have to recheck here.
    }
    unless (exists $$test{'SNPS'}{$rsid}){
        push @missing, $rsid;
    }
}

if (scalar @missing > 0) {
    print "The following " . scalar @missing . " SNPs have not been analysed because they were not found in the 1000 genomes phase 1 integrated call data\n";
    print join("\n", @missing) . "\n";
}
if (defined $ld) {
    if ($output < $input) {
        say "For $label, $input SNPs provided, " . scalar @snps . " retained, " . scalar @missing . " not analysed, "  . scalar(keys %$ld_excluded) . " LD filtered at $ld,";
    }
}

# only pick background snps matching snps that had bitstrings originally.
my @foundsnps = keys %{$$test{'SNPS'}};
my $snpcount = scalar @foundsnps;
print "Test SNPs analysed $snpcount\n";

# identify the gc, maf and tss, and then make bkgrd picks
my $picks = match(\%$test, $bkgd, $datadir, $per, $reps) unless defined $overlap;

# for bkgrd set need to get distribution of counts instead
# make a hash of data -> cell -> bkgrd-Set -> overlap counts
my %bkgrd; #this hash is going to store the bkgrd overlaps

# Get the bits for the background sets and process
my $backsnps;

if (defined $picks){
    foreach my $bkgrd (keys %{$picks}){
        #$rows = get_bits(\@{$$picks{$bkgrd}}, $sth);
        $rows = get_bits(\@{$$picks{$bkgrd}}, $dbh);
        $backsnps += scalar @$rows; #$backsnps is the total number of background SNPs analysed
        unless (scalar @$rows == scalar @foundsnps){
                print "Background " . $bkgrd . " only " . scalar @$rows . " SNPs out of " . scalar @foundsnps . "\n";
            }
            my $result = process_bits($rows, $cells, $data);
            foreach my $cell (keys %{$$result{'CELLS'}}){
                push @{$bkgrd{$cell}}, $$result{'CELLS'}{$cell}{'COUNT'}; # accumulate the overlap counts by cell
            }
            if (defined $bkgrdstat){
                bkgrdstat($result, $lab, $bkgrd);
            }
    }
}

$dbh->disconnect();

#Having got the test overlaps and the bkgd overlaps now calculate Zscores, pvalue and output the table to be read into R for plotting.
my $time = time(); # time is used to label the output directories.
my $resultsdir;
if (defined $peaks){
    $resultsdir = "$cwd/$lab.peaks.$time";
}
else{
    $resultsdir = "$cwd/$lab.$time";
}
mkdir $resultsdir;
my $filename = "$lab.chart.tsv";
open my $ofh, ">", "$resultsdir/$filename" or die "Cannot open $resultsdir/$filename: $!"; #should grab a process number for unique name here
if (defined $overlap){
    print $ofh join("\t", "Cell", "Tissue", "File", "SNPs", "Number", "Accession") ."\n";
}
else{
    print $ofh join("\t", "Zscore", "Pvalue", "Cell", "Tissue", "File", "SNPs", "Number", "Accession") ."\n";
}
my $n = 1;
my $pos = 0;

my %tissuecount;
foreach my $cell (keys %$tissues){
    my $tissue = $$tissues{$cell}{'tissue'};
    $tissuecount{$tissue}++;
}

my $tissuecount = scalar keys %tissuecount;

$t1 = $t1/$tissuecount; # bonferroni correction by number of tissues
$t2 = $t2/$tissuecount;

$t1 = -log10($t1);
$t2 = -log10($t2);

open my $bfh, ">", "background.tsv" or die "Cannot open background.tsv";

foreach my $cell (sort {ncmp($$tissues{$a}{'tissue'},$$tissues{$b}{'tissue'}) || ncmp($a,$b)} @$cells){ # sort by the tissues alphabetically (from $tissues hash values)
    # ultimately want a data frame of names(results)<-c("Zscore", "Cell", "Tissue", "File", "SNPs")

    my $snp_string = "";
    $snp_string = join(",", @{$$test{'CELLS'}{$cell}{'SNPS'}}) if defined $$test{'CELLS'}{$cell}{'SNPS'}; # This gives the list of overlapping SNPs for use in the tooltips. If there are a lot of them this can be a little useless
    my ($shortcell, undef) = split('\|', $cell); # undo the concatenation from earlier to deal with identical cell names.

    my $teststat = ($$test{'CELLS'}{$cell}{'COUNT'} or 0); #number of overlaps for the test SNPs

    # binomial pvalue, probability of success is derived from the background overlaps over the tests for this cell
    # $backsnps is the total number of background SNPs analysed
    # $tests is the number of overlaps found over all the background tests
    if (defined $overlap){
        print $ofh join("\t", $shortcell, $$tissues{$cell}{'tissue'}, $$tissues{$cell}{'file'}, $snp_string, $n, $$tissues{$cell}{'acc'}) . "\n";
    }
    else{
        say $bfh join("\t", @{$bkgrd{$cell}});
        my $tests;
        foreach (@{$bkgrd{$cell}}){
            $tests+= $_;
        }
        my $p = sprintf("%.6f", $tests/$backsnps);

        # binomial probability for $teststat or more hits out of $snpcount snps
        # sum the binomial for each k out of n above $teststat
        my $pbinom;
        foreach my $k ($teststat .. $snpcount){
            $pbinom += binomial($k, $snpcount, $p);
        }
        if ($pbinom >1) {
            $pbinom = 1;
        }
        $pbinom = -log10($pbinom);
        # Z score calculation
        my $mean = mean(@{$bkgrd{$cell}});
        my $sd = std(@{$bkgrd{$cell}});
        my $zscore;
        if ($sd == 0){
            $zscore = "NA";
        }
        else{
            $zscore = sprintf("%.3f", ($teststat-$mean)/$sd);
        }
        if ($pbinom >=$t2){
            $pos++;
        }
        print $ofh join("\t", $zscore, $pbinom, $shortcell, $$tissues{$cell}{'tissue'}, $$tissues{$cell}{'file'}, $snp_string, $n, $$tissues{$cell}{'acc'}) . "\n";
    }
        $n++;
}


# fdr calculation isn't valid currently
#my $fdr = fdr($pos, $snpcount, $cellcount);
#say "$filename\t$pos positive lines at FDR = $fdr at Z >= 3.39";

unless (defined $noplot){
    #Plotting and table routines
    unless (defined $overlap){
        Chart($filename, $lab, $resultsdir, $tissues, $cells, $label, $t1, $t2, $data); # basic pdf plot
        dChart($filename, $lab, $resultsdir, $data, $label, $t1, $t2); # rCharts Dimple chart
    }
    table($filename, $lab, $resultsdir, $overlap); # Datatables chart
}
