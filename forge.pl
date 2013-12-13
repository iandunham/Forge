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

In each of the graphics the colouring should be consistent. Blue (Z < 2.58), light red or pink (2.58 =< Z < 3.39), red or dark red (Z >= 3.39 ) for the 99% and 99.9% cIs. Or whatever other thresholds are specified.

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

Specify the minimum number of SNPs to be allowed. Default is 20.

=item B<thresh>

Alter the default Z score thresholds. Give a comma separate list of two e.g.2.58,3.39 for the defaults

=item B<format>

if f is specified, specify the file format as follow:

rsid = list of snps as rsids each on a separate line

bed  = File given is a bed file of locations (chr\tbeg\tend) aka Personal Genome SNP format.  bed format should be 0 based and the chromosome should be given as chrN. Hoever will also accept chomosomes as just N (ensembl) and 1-based format where beg and end are the same

vcf = File given is a vcf file.

tabix = File contains SNPs in tabix format.

ian = 1-based chr\tbeg\tend\trsid\tpval

=item B<bkgrd>

Output background stats for investigation.

=item B<reps>

THe number of background matching sets to pick and analyse. Default 100.

=item B<noplot>

Just make the data file, don't plot.

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
use Storable;
use Getopt::Long;
use File::Basename;
use Config::IniFiles;
use Pod::Usage;
use Scalar::Util qw(looks_like_number);

my $cwd = getcwd;

my ($bkgd, $data, $peaks, $label, $file, $format, $min_snps, $bkgrdstat, $noplot, $reps, $help, $man, $thresh, @snplist);

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
    'help|h|?'   => \$help,
    'man|m'      => \$man,

);

pod2usage(1) if ($help);
pod2usage(-verbose => 2) if ($man);

unless (defined $min_snps){
    $min_snps = 20; # the minimum number of snps allowed for test.
}

unless (defined $data ){
    $data = "erc"; # define which data we are dealing with for the bitstrings.
}
unless (defined $label){
    $label = "No label given";
}
(my $lab = $label) =~ s/\s/_/g;
$lab = "$lab.$data";

unless (defined $format){
    $format = 'rsid';
}
my $dirname = dirname(__FILE__);

my $cfg = Config::IniFiles->new( -file => "$dirname/forge.ini" );
my $datadir = $cfg->val('Files', 'datadir');

# percentile bins for the bkgrd calculations. This is hard coded so there are enough SNPs to choose from, but could later be altered.
my $per = 10;
# number of sets to analyse for bkgrd. Again currently hardcoded to 100

unless (defined $reps){
    $reps = 100;
}

unless (defined $bkgd){
    $bkgd = "gwas";
}
my ($t1, $t2);

if (defined $thresh){
    ($t1, $t2) = split(",", $thresh);
    unless (looks_like_number($t1) && looks_like_number($t2)){
        die "You must specify numerical thersholds in a comma separated list";
    }
}
else{
    $t1 = 2.58;
    $t2 = 3.39;
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
    # would be better if this was module or sub :-(
    my $sth = $dbh->prepare("SELECT rsid FROM bits WHERE location = ?");
    open my $fh, "<", $file or die "cannot open file $file : $!";
    if ($format =~ /rsid/){
        while (<$fh>){
            chomp;
            my @rsid = split /\:/, $_;
            my $rsid = pop @rsid; # take the last one for want of a better idea.
            push @snps, $rsid;
        }
    }
    elsif ($format =~ /ian/){
        while (<$fh>){
            my ($chr, $beg, $end, $rsid, undef) = split "\t", $_;
            my @rsid = split /\:/, $rsid; # to deal with multiple rsids
            $rsid = pop @rsid; # take the last one for want of a better idea. Can't take all as they are the same thing.
            push @snps, $rsid;
        }
    }
    elsif ($format =~ /vcf/){
        while (<$fh>){
            next if /^#/;
            my ($chr, $beg, $rsid) = split "\t", $_;
            unless ($chr =~ /^chr/){
                $chr = "chr". $chr;
            }
            if ($rsid =~/^rs\d+/){
                push @snps, $rsid;
            }
            else {
                my $loc = "$chr:$beg-$beg";
                #get the rsid from the db
                $rsid = fetch_rsid($loc, $sth);
                push @snps, $rsid if defined $rsid;
            }
        }
    }
    elsif ($format =~ /bed|tabix/){
        while (<$fh>){
            my $loc;
            if ($format =~/bed/){
                next if /^track/;
                my ($chr, $beg, $end) = split "\t", $_;
                unless ($chr =~ /^chr/){
                    $chr = "chr". $chr;
                }
                $loc = "$chr:$end-$end";
            }
            elsif ($format =~ /tabix/){
                chomp;
                $loc = $_;
            }
            #get the $rsid from the db
            my $rsid = fetch_rsid($loc, $sth);
            push @snps, $rsid if defined $rsid;
        }
    }
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


# Check we have enough SNPs
if (scalar @snps < $min_snps){
    pod2usage(-verbose => 2, -message => "Fewer than $min_snps SNPs provided. Analysis not run\n\n", -noperldoc => 1);
}

# Connect to the sqlite database file which contains the tables for each data
#my $sth = $dbh->prepare("SELECT * FROM bits WHERE rsid IN (?)"); # ideally we could prepare once with the right number of placeholders and then execute on each SNP set, providing they were the same length. Problem is that sometimes fail to find a match for 1 SNP.

# get the cell list array and the hash that connects the cells and tissues
my ($cells, $tissues) = get_cells($data, $dbh);

# get the bit strings for the test snps from the database file
#my $rows = get_bits(\@snps, $sth);
my $rows = get_bits(\@snps, $dbh);

# unpack the bitstrings and store the overlaps by cell.
my $test = process_bits($rows, $cells, $data);
if (defined $bkgrdstat){
    open my $bfh, ">", "$lab.bkgrd.stats" or die "cannot open $lab.bkgrd.stats";
        my (@maf, @tss, @gc);
        foreach my $rsid (keys %{$$test{'SNPS'}}){
            my ($maf, $tss, $gc) = split "\t", $$test{'SNPS'}{$rsid}{'PARAMS'};
            push @maf, $maf;
            push @tss, $tss;
            push @gc, $gc;
        }
    say $bfh join("\t", "test", "maf", @maf);
    say $bfh join("\t", "test", "tss", @tss);
    say $bfh join("\t", "test", "gc", @gc);
}


# Identify SNPs that weren't found and warn about them.
my @missing;
foreach my $rsid (@snps){
    unless (exists $$test{'SNPS'}{$rsid}){
        push @missing, $rsid;
    }
}
if (scalar @missing > 0) {
    print "The following " . scalar @missing . " SNPs have not been analysed because they were not found in the 1000 genomes phase 1 integrated call data\n";
    print join("\n", @missing) . "\n";
}

# only pick background snps matching snps that had bitstrings originally.
my @foundsnps = keys %{$$test{'SNPS'}};
my $snpcount = scalar @foundsnps;
print "Test SNPs analysed $snpcount\n";

# identify the gc, maf and tss, and then make bkgrd picks
my $picks = match(\%$test, $bkgd);

# for bgrd set need to get distribution of counts instead
# make a hash of data -> cell -> bkgrd-Set -> overlap counts
my %bkgrd; #this hash is going to store the bkgrd overlaps

# Get the bits for the background sets and process

foreach my $bkgrd (keys %{$picks}){
    #$rows = get_bits(\@{$$picks{$bkgrd}}, $sth);
    $rows = get_bits(\@{$$picks{$bkgrd}}, $dbh);
    unless (scalar @$rows == scalar @foundsnps){
        print "Background " . $bkgrd . " only " . scalar @$rows . " SNPs out of " . scalar @foundsnps . "\n";
    }
    my $result = process_bits($rows, $cells, $data);
    foreach my $cell (keys %{$$result{'CELLS'}}){
        push @{$bkgrd{$cell}}, $$result{'CELLS'}{$cell}{'COUNT'}; # accumulate the overlap counts by cell
    }
    if (defined $bkgrdstat){
        open my $bfh, ">>", "$lab.bkgrd.stats" or die "cannot open $lab.bkgrd.stats";
        my (@maf, @tss, @gc);
        foreach my $rsid (keys %{$$result{'SNPS'}}){
            my ($maf, $tss, $gc) = split "\t", $$result{'SNPS'}{$rsid}{'PARAMS'};
            push @maf, $maf;
            push @tss, $tss;
            push @gc, $gc;
        }
        say $bfh join("\t", $bkgrd, "maf", @maf);
        say $bfh join("\t", $bkgrd, "tss", @tss);
        say $bfh join("\t", $bkgrd, "gc", @gc);
    }
}

$dbh->disconnect();

#Having got the test overlaps and the bkgd overlaps now calculate Zscores and output the table to be read into R for plotting.
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
print $ofh join("\t", "Zscore", "Cell", "Tissue", "File", "SNPs", "Number", "Accession") ."\n";
my $n =1;

my $pos = 0;

open my $bfh, ">", "background.tsv" or die "Cannot open background.tsv";

foreach my $cell (sort {ncmp($$tissues{$a}{'tissue'},$$tissues{$b}{'tissue'}) || ncmp($a,$b)} @$cells){ # sort by the tissues alphabetically (from $tissues hash values)
    # ultimately want a data frame of names(results)<-c("Zscore", "Cell", "Tissue", "File", "SNPs")
    say $bfh join("\t", @{$bkgrd{$cell}});
    my $mean = mean(@{$bkgrd{$cell}});
    my $sd = std(@{$bkgrd{$cell}});
    my $teststat = $$test{'CELLS'}{$cell}{'COUNT'};
    my $zscore;
    if ($sd == 0){
        $zscore = "NA";
    }
    else{
        $zscore = sprintf("%.3f", ($teststat-$mean)/$sd);
    }
    if ($zscore >=$t2){
        $pos++;
    }
    my $snp_string = "";
    $snp_string = join(",", @{$$test{'CELLS'}{$cell}{'SNPS'}}) if defined $$test{'CELLS'}{$cell}{'SNPS'}; # This gives the list of overlapping SNPs for use in the tooltips. If there are a lot of them this can be a little useless
    my ($shortcell, undef) = split('\|', $cell); # undo the concatenation from earlier to deal with identical cell names.
    print $ofh join("\t", $zscore, $shortcell, $$tissues{$cell}{'tissue'}, $$tissues{$cell}{'file'}, $snp_string, $n, $$tissues{$cell}{'acc'}) . "\n";
    $n++;
}

my $cellcount = scalar @$cells;
my $fdr = fdr($pos, $snpcount, $cellcount);
say "$filename\t$pos positive lines at FDR = $fdr";

unless (defined $noplot){
    #Plotting and table routines
    Chart($filename, $lab, $resultsdir); # basic pdf plot
    #rChart($filename, $lab, $resultsdir); # rCharts polychart plot
    dChart($filename, $lab, $resultsdir); # rCharts Dimple chart
    table($filename, $lab, $resultsdir); # Datatables chart
    #hChart("$lab.chart.tsv", $label);
}

### Subroutines ###

sub match{
    # identifies the bins that each of the snps lies in, and then picks 100 matching SNPs.
    my $snps = shift; # ahash that contians the found snps
    my $bkgd = shift;
    my ($bins, $params, %bins, %params);
    # load up the stored hashes that contain the bins of snps by gc, maf, and tss distance. There is one for each of the bkgd set possibilities(gwas or omni).
    # These are precalculated according to the parameters that are hard coded above.
    # the hash to load is defined by the bkgd option - defualts to 'gwas'
    if ($bkgd =~ /omni/i){
        $bins = $datadir . "omni.snp_bins.$per";
        $params = $datadir . "omni.snp_params.$per";
    }
    else{
        $bins = $datadir . "snp_bins.$per";
        $params = $datadir . "snp_params.$per";
    }
    if (-e $bins && -e $params){
        %bins = %{ retrieve($bins) };
        %params = %{ retrieve($params)};
    }
    my (%picks);

    foreach my $rs (keys %{$$snps{'SNPS'}}){
        srand;
        my ($maf, $tss, $gc) = split("\t", join("\t", $$snps{'SNPS'}{$rs}{'PARAMS'}));
        #$rs is the test snp, $rsid is the matched snp.
        my ($i, $j, $k) = assign ($gc, $tss, $maf, \%params);

        my $range = scalar @{$bins{$i}{$j}{$k}};
        for (my $n = 1; $n <= $reps; $n++) {
            my ($snp_string, $rsid);
            while (1){
                my $pick = int(rand($range));
                $snp_string = ${$bins{$i}{$j}{$k}}[$pick]; #pick the $pick'th element in the array as the chosen SNP "
                (undef, undef, undef, $rsid) = split /\t/,  $snp_string;
                last unless $rsid eq $rs; # must not pick the test snp itself.
            }
            push @{$picks{$n}}, $rsid; # each $n array is a set of snps matching the test set/ it is allowed to pick the same SNP more than once in this backgrouns selection
        }
    }
    return \%picks;
}

sub process_bits{
    # Processes the bitstrings to get a count of overlaps for each cell type.
    my ($rows, $cells, $data) = @_;
    my %test;
    foreach my $row (@{$rows}){
        my ($location, $rsid, $sum, $bit, $maf, $tss, $gc);
        if ($data eq "erc"){
            ($location, $rsid, undef, undef, $bit, $sum, $maf, $tss, $gc) =  @$row;
        }
        else{
            ($location, $rsid, $bit, $sum, undef, undef, $maf, $tss, $gc) = @$row;
        }
        $test{'SNPS'}{$rsid}{'SUM'} = $sum;
        $test{'SNPS'}{$rsid}{'PARAMS'} = join("\t", $maf, $tss, $gc);
        my @bits = split "", $bit;
        my $index = 0;
        foreach my $cell (@$cells){
            $test{'CELLS'}{$cell}{'COUNT'} += $bits[$index];
            push @{$test{'CELLS'}{$cell}{'SNPS'}}, $rsid if $bits[$index] == 1;
            $index++;
        }
    }
    return \%test;
}

sub get_bits{
    #get the bitstrings for an array of snps from the sqlite db
    my ($snps, $dbh) = @_;
    my @results;
    my $args = join ("','", @$snps);
    my $sth = $dbh->prepare("SELECT * FROM bits WHERE rsid IN ('$args')");
    $sth->execute();
    my $result = $sth->fetchall_arrayref();
    $sth->finish();
    foreach my $row (@{$result}){
      push @results, $row;
    }
    return \@results;# return the bitstring line from the database
}


sub fetch_rsid{
    #gets the rsid for a SNP where a location is given
    my ($loc, $sth) = @_;
    $sth->execute($loc);
    my $result = $sth->fetchall_arrayref();
    my $rsid;
    foreach my $row (@{$result}){
        $rsid = $$row[0];
    }
    $sth->finish();
    if (defined $rsid &&$rsid =~ /^rs\d+/){
        return $rsid;
    }
    else{
        return "no RSID match for $loc";
    }
}

sub get_cells{
    # read the correct cell list based on data (erc -encode). Also gets the tissue names for the cells.
    my $data = shift;
    my $dbh = shift;
    my $table = join('_', "cells", $data);
    my $sth = $dbh->prepare("SELECT shortcell,tissue,file,acc FROM $table");
    $sth->execute();
    my $ver = $sth->fetchall_arrayref();
    $sth->finish();
    my ($cells, $tissues, $acc);
    foreach my $row (@$ver){
        my $cell = shift @$row;
        my $tissue = shift @$row;
        my $file = shift @$row;
        my $acc = shift @$row;
        $cell = "$cell|$file"; # Sometimes the same cell is used twice, with a differnt file so need to record separately (e.g. WI-38).
        push @$cells, $cell;
        $$tissues{$cell}{'tissue'} = $tissue; # this is the hash that is used to connect cells and tissues and ultimately provide the sorting order
        $$tissues{$cell}{'file'} = $file;
        $$tissues{$cell}{'acc'} = $acc;
    }
    #print Dumper %$tissues;
    return ($cells, $tissues); # return
}

sub assign{
    #sub routine to assign any maf, gc, tss values to the percentile bins
    my ($gc, $tss, $maf, $params) = @_;
    my ($i, $j, $k);
    my $n = 1;
    foreach my $pc (@{$$params{'gc'}}){
        if ($gc <= $pc) {
            $i = $n;
        }
        else{
            $n++;
        }
    }
    $n=1;
    foreach my $pc (@{$$params{'tss'}}){
        if ($tss <= $pc) {
            $j = $n;
        }
        else{
            $n++;
        }
    }
    $n=1;
    foreach my $pc (@{$$params{'maf'}}){
        if ($maf <= $pc) {
            $k = $n;
        }
        else{
            $n++;
        }
    }
    return ($i, $j, $k);
}

sub mean {
    # calculates the biased mean of an array
    #
    # pass it a float array and it will return the mean
    # reused from Ben Brown
    my $sum = 0;
    foreach (@_){
        $sum+= $_;
        }
    return $sum/($#_+1);
}

sub var {
    # calculates the biased variance of an array
    #
    # pass it a float array and it will return the variance
    # reused from Ben Brown
    my $ev = mean(@_);
    my $sum = 0;
    foreach (@_) { $sum += ($_ - $ev)**2 };

    return $sum/($#_+1);
}

# calulates the standard deviation of an array: this is just the sqrt of the var
sub std { sqrt(var(@_)) }

sub fdr{
    my ($tp, $snps, $cells) = @_;
    if ($tp == 0){
        return "NA";
    }
    else{

        my $fpr = 0.0085 * exp(-0.04201 * $snps) + 0.00187; # from simulations of random data  0.0085*exp(-0.04201. SNPs) + 0.00187
        my $fdr = ($cells * $fpr) / $tp;
        return $fdr;
    }
}

sub Chart{
    # This is the original code using standard R plot to generate a static pdf.
    print "Making static chart.\n";
    my ($filename, $lab, $resultsdir) = @_;
    my $Rdir = $resultsdir;
    my $chart = "$lab.chart.pdf";
    my $rfile = "$Rdir/$lab.chart.R";

    # make plot, first calculate where dividing lines are:
    my (@lines, @label_pos, @labels, @tissue_txt);
    my $n =1;
    my $last = '0';
    foreach my $cell (sort {ncmp($$tissues{$a}{'tissue'},$$tissues{$b}{'tissue'}) || ncmp($a,$b)} @$cells) {
        my $tissue = $$tissues{$cell}{'tissue'};
        push @tissue_txt, $tissue;
        unless ($tissue eq $last){
            push @lines, $n-0.5;
            push @labels, $tissue;
        }
        $n++;
        $last = $tissue;
    }

    my $length = scalar @tissue_txt;
    my $index =0;
    foreach my $value (@lines) {
        if (defined $lines[$index+1]){
            my $pos = $value + ($lines[$index+1]-$value)/2;
            push @label_pos, $pos;
        }
        else {
            my $pos = $value + ($length-$value)/2;
            push @label_pos, $pos;
        }
        $index++;
    }

    open my $rfh, ">", "$rfile";
    print $rfh "setwd(\"$Rdir\")
results<-read.table(\"$filename\",header=TRUE,sep=\"\t\")
results\$Class<-cut(results\$Zscore, breaks =c(min(results\$Zscore), $t1, $t2, max(results\$Zscore)), labels=FALSE, include.lowest=TRUE) # 99.9 and 99% CIs 1, 2, 3
pdf(\"$chart\", width=22.4, height=7)
palette(c(\"steelblue3\",\"pink2\",\"red\"))
ymin1 = min(results\$Zscore, na.rm=TRUE)*1.1
ymax1 = max(results\$Zscore, na.rm=TRUE)*1.1
ymax = max(c(abs(ymin1),ymax1))
ymin = -ymax
par(mar=c(9,4,3,1)+0.1)
plot(results\$Zscore,ylab=\"Z score\",xlab=\"\",main=\"Proportion of SNPs, DNase1 sites (probably TF sites) which are present in cell lines for $label\",ylim=c(ymin,ymax), las=2, las=2, pch=19,col=results\$Class, xaxt='n')
axis(1, seq(1,length(results\$Cell)),labels=results\$Cell, las=2, cex.axis=0.7)
mtext(1,text=\"Cell\",line=7,cex=1.2)
#abline(h=-$t1, col=\"lightpink1\") # Z score of 2.58 = 99 % probability
abline(h=$t1, col=\"lightpink1\")
#abline(h=-$t2, col=\"lightpink1\", lty=2)
abline(h=$t2, col=\"lightpink1\", lty=2)
text(c(-1),$t1+0.2,c(\"Z = $t1\"),col=\"lightpink1\",adj=1,cex=0.8)
#text(c(-1),-$t1+0.16,c(\"1%\"),col=\"lightpink1\",adj=1,cex=0.8)
text(c(-1),$t2+0.2,c(\"Z = $t2\"),col=\"lightpink1\",adj=1,cex=0.8)
#text(c(-1),-$t2+0.16,c(\"0.1%\"),col=\"lightpink1\",adj=1,cex=0.8)
palette(\"default\")\n";

    foreach my $pos (@lines){
        print $rfh "lines(c($pos,$pos),c(-22,22),lty=6,col=\"#00000070\")\n" unless $pos == 0.5;
    }
    $index = 0;
    foreach my $tissue (@labels){
        print $rfh "text(c(" . $label_pos[$index] . "),ymax,c(\"" . $tissue . "\"),col=\"burlywood3\",adj=1,srt=90,cex=0.8)\n";
        $index++;
    }
    print $rfh "dev.off()\n";
#run the R code
    system "R --no-save --quiet --slave < $rfile";
}

sub rChart{
    # Makes a polcharts plot : note X axis labelling is problematical
    print "Making rChart.\n";
    my ($filename, $lab, $resultsdir) = @_;
    my $chart = "$lab.rchart.html";
    my $Rdir = $resultsdir;
    my $rfile = "$Rdir/$lab.rChart.R";
    open my $rcfh, ">", "$rfile";
    print $rcfh "setwd(\"$Rdir\")
results<-read.table(\"$filename\", header = TRUE, sep=\"\t\")
results\$Colour<- 0 + (results\$Zscore < $t2) + (results\$Zscore < $t1)  # 99.9 and 99% CIs
require(rCharts)
r1 <- rPlot(Zscore ~ Cell, data=results, color=\"bin(Colour, 0.25)\", type=\"point\", tooltip = \"function(item){ return (item.Zscore + '\\\\n' + item.Cell + '\\\\n' + item.Tissue + '\\\\n' + item.File + '\\\\n' + item.SNPs + '\\\\n' + item.Accession + '\\\\n')}\")
#r1\$guides(color=list(scale = list(type = \'gradient\', lower = \'\#CCC\', upper = \'\#000\'))) # optional code to make a grey scale
r1\$addParams(width = 2000, height=600, title=\"$label overlaps with $data DHS\")
ymin1 = min(results\$Zscore, na.rm=TRUE)*1.2
ymax1 = max(results\$Zscore, na.rm=TRUE)*1.2
ymax = max(c(abs(ymin1),ymax1))
ymin = -ymax
r1\$guides(x = list(numticks = length(unique(results\$Cell)), levels=results\$Cell), y = list(min = ymin, max = ymax))
r1\$save('$chart', cdn = F)
##r1\$show() #makes a temp file\n";

system "R --no-save --quiet --slave < $rfile";
}

sub dChart{
    # Make dimple interactive chart.
    print "Making dChart.\n";
    my ($filename, $lab, $resultsdir) = @_;
    my $chart = "$lab.dchart.html";
    my $Rdir = $resultsdir;
    my $rfile = "$Rdir/$lab.dChart.R";
   open my $rcfh, ">", "$rfile";
    print $rcfh "setwd(\"$Rdir\")
results<-read.table(\"$filename\", header = TRUE, sep=\"\t\")
results\$Class<-cut(results\$Zscore, breaks =c(min(results\$Zscore), $t1, $t2, max(results\$Zscore)), labels=FALSE, include.lowest=TRUE) # 99.9 and 99% CIs 1, 2, 3
require(rCharts)
d1 <- dPlot(
  y = \"Zscore\",
  x = c(\"Cell\", \"Tissue\", \"SNPs\", \"Number\", \"Accession\"),
  groups = \"Class\",
  data = results,
  type = \"bubble\",
  width = 2000,
  height = 1500,
  bounds = list(x=90,y=50,height=600,width=1850),
  id = \"chart.$lab\"
)\n";
if ($data =~ /erc/){
    print $rcfh "d1\$xAxis( type = \"addCategoryAxis\", grouporderRule = \"Tissue\", orderRule = \"Cell\")\n";
}
else {
    print $rcfh "d1\$xAxis( type = \"addCategoryAxis\", grouporderRule = \"Tissue\", orderRule = \"Number\")\n";
}

print $rcfh "d1\$yAxis( type = \"addMeasureAxis\" )
d1\$colorAxis(
   type = \"addColorAxis\",
   colorSeries = \"Class\",
   palette = c(\"lightblue\",\"pink\",\"red\") )
d1\$addParams(title=\"$label overlaps with $data DHS\")
d1\$save('$chart', cdn = F)\n";

system "R --no-save --quiet --slave < $rfile";
}

sub table{
    # Make Datatables table
    print "Making Table.\n";
    my ($filename, $lab, $resultsdir) = @_;
    my $chart = "$lab.table.html";
    my $Rdir = $resultsdir;
    my $rfile = "$Rdir/$lab.table.R";
    open my $rcfh, ">", "$rfile";
    print $rcfh "setwd(\"$Rdir\")
    data<-read.table(\"$filename\", header = TRUE, sep=\"\t\")
    results<-data.frame(data\$Cell, data\$Tissue, data\$Accession, data\$Zscore, data\$SNPs)
    names(results)<-c(\"Cell\", \"Tissue\", \"Accession\", \"Zscore\", \"SNPs\")
    require(rCharts)
    dt <- dTable(
      results,
      sScrollY= \"600\",
      bPaginate= F,
      sScrollX= \"100%\",
      sScrollXInner= \"110%\"
    )
    dt\$save('$chart', cdn = F)";
    system "R --no-save --quiet --slave < $rfile";
}

# highcharts interface has problems with plotting
#sub hChart{
#    print "Making hChart.\n";
#    my ($filename, $label) = @_;
#    my $chart = $lab . ".hchart.html";
#    my $Rdir = $cwd;
#    my $rfile = $lab . ".hChart.R";
#    open my $rcfh, ">", "$Rdir/$rfile";
#    print $rcfh "setwd(\"$Rdir\")
#results<-read.table(\"$filename\", header = TRUE, sep=\"\t\")
#results\$Class<-cut(results\$Zscore, breaks =c(min(results\$Zscore), $t1, $t2, max(results\$Zscore)), labels=FALSE, include.lowest=TRUE) # 99.9 and 99% CIs 1, 2, 3
#require(rCharts)
#h1 <- hPlot(Zscore ~ Number, data=results, type=\"scatter\", radius=5)
#h1\$addParams(width = 2000, height=600, title=list(\"$label overlaps with $data DHS\"))
#h1\$xAxis(title = list(text = \"Cell\"), labels = list(rotation=-90, align=\"right\"), categories = results\$Cell\)
#h1\$save('$chart', cdn = F)";
#
#system "R --no-save --quiet --slave < $rfile";
#}
