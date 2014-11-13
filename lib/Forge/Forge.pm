package Forge::Forge;

use 5.010;
use strict;
use warnings FATAL => 'all';
use Storable;

=head1 NAME

Forge - The great new Forge!

=head1 VERSION

Version 0.01

=cut

my $MAX_SQL_VARIABLES = 999;
our $VERSION = '0.01';
our (@ISA, @EXPORT);
use Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(process_file match process_bits get_bits fetch_rsid fetch_loc get_cells assign bkgrdstat ld_filter);

=head1 SYNOPSIS

Provide functional modules for Forge

=head1 EXPORT

process_file
match
process_bits
get_bits
fetch_rsid
fetch_loc
get_cells
assign
bkgrdstat

=head1 SUBROUTINES/METHODS

=head2 process_file

Processes various file formats.

=cut

sub bkgrdstat{
    my ($test, $lab, $flag) = @_;
    my $bfh;
    my $file = "$lab.bkgrd.stats";
    if ($flag eq "test") {
        open $bfh, ">", $file or die "cannot open $file";
    }
    else{
        open $bfh, ">>", $file or die "cannot open $file";
    }
    my (@maf, @tss, @gc);
    foreach my $rsid (keys %{$$test{'SNPS'}}){
        my ($maf, $tss, $gc) = split "\t", $$test{'SNPS'}{$rsid}{'PARAMS'};
        push @maf, $maf;
        push @tss, $tss;
        push @gc, $gc;
    }
    say $bfh join("\t", $flag, "maf", @maf);
    say $bfh join("\t", $flag, "tss", @tss);
    say $bfh join("\t", $flag, "gc", @gc);
}

=head2 process_file

Processes various file formats.

=cut

sub process_file {
    my ($fh, $format, $sth, $filter) = @_;
    my @snps;
    if ($format =~ /rsid/){
        while (<$fh>){
            chomp;
            my $rs;
            if (defined $filter) {
                my $pval;
                ($rs, $pval) = split /\s+/, $_;
                next unless $pval >= $filter;
            }
            else{
                ($rs, undef) = split /\s+/, $_; # remove anything that is not supposed to be there :-)
            }
            my @rsid = split /\:/, $rs;
            my $rsid = pop @rsid; # take the last one for want of a better idea.
            push @snps, $rsid;
        }
    }
    elsif ($format =~ /ian/){
        while (<$fh>){
            my ($chr, $beg, $end, $rsid, $p, $pval) = split "\t", $_;
            if (defined $filter) {
                next unless $pval >= $filter;
            }
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
            chomp;
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
    return @snps;
}

=head2 match

Identifies the bins that each of the snps in a snp hash lies in, and then picks matching SNPs for the number of reps specified.

=cut

sub match{
    my ($snps, $bkgd, $datadir, $per, $reps) = @_;
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

=head2 process_bits

Processes bitstrings to get a count of overlaps for each cell type.

=cut

sub process_bits{
    my ($rows, $cells, $data) = @_;
    my %test;
    foreach my $row (@{$rows}){
        my ($location, $rsid, $sum, $bit_string, $maf, $tss, $gc);
        if ($data eq "erc"){
            ($location, $rsid, undef, undef, $bit_string, $sum, $maf, $tss, $gc) =  @$row;
        }
        else{
            ($location, $rsid, $bit_string, $sum, undef, undef, $maf, $tss, $gc) = @$row;
        }
        $test{'SNPS'}{$rsid}{'SUM'} = $sum;
        $test{'SNPS'}{$rsid}{'PARAMS'} = join("\t", $maf, $tss, $gc);
        die if (scalar(@$cells) ne length($bit_string));
        my $index = 0;
        foreach my $cell (@$cells){
            ## $bit_string is a string made of 0s and 1s. If it is a 1 for this position, count and push
            if (substr($bit_string, $index, 1)) {
                $test{'CELLS'}{$cell}{'COUNT'}++;
                push @{$test{'CELLS'}{$cell}{'SNPS'}}, $rsid;
            }
            $index++;
        }
    }
    return \%test;
}

=head2 get_bits

Get the bitstrings for an array of snps from the sqlite db

=cut

sub get_bits{

    my ($snps, $dbh) = @_;
    my @results;


    for (my $loop = 0; $loop * $MAX_SQL_VARIABLES < @$snps; $loop++) {
        my $start = $loop * $MAX_SQL_VARIABLES;
        my $end = ($loop + 1) * $MAX_SQL_VARIABLES - 1;
        $end = @$snps -1 if ($end >= @$snps);

        my $sql = "SELECT * FROM bits WHERE rsid IN (?". (",?" x ($end - $start)).")";
        my $sth = $dbh->prepare($sql); #get the blocks form the ld table
        $sth->execute(@$snps[$start..$end]);

        my $result = $sth->fetchall_arrayref();
        $sth->finish();
        foreach my $row (@{$result}){
            push @results, $row;
        }
    }

    return \@results;# return the bitstring line from the database
}

=head2 fetch_rsid

gets the rsid for a SNP where a location is given.

=cut

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

=head2 fetch_loc

Filter SNPs from the SNP list if they are in LD at r2>= 0.8 or r2 >=0.1.  THe rationale is that the first SNP to be indentified in a block is chosen, and others are removed.

=cut

sub ld_filter{
    #@snps = ld_filter(@snps, $ld, $dbh);
    my ($snps, $r2, $dbh) = @_;
    my %ld_excluded; # a hash to store SNPs found in LD with a SNP in the list
    my @snps_filtered; # The list of SNPs filtered
    my %snps;
    foreach my $snp (@$snps){
        $snps{$snp} = 1;
    }
    for (my $loop = 0; $loop * $MAX_SQL_VARIABLES < @$snps; $loop++) {
        my $start = $loop * $MAX_SQL_VARIABLES;
        my $end = ($loop + 1) * $MAX_SQL_VARIABLES - 1;
        $end = @$snps -1 if ($end >= @$snps);

        my $sql = "SELECT rsid,$r2 FROM ld WHERE rsid IN (?". (",?" x ($end-$start)).")";
        my $sth = $dbh->prepare($sql); #get the blocks form the ld table
        $sth->execute(@$snps[$start..$end]);
        my $result = $sth->fetchall_arrayref();
        $sth->finish();
        foreach my $row (@{$result}){
            my ($snp, $block) = @$row;
            next if exists $ld_excluded{$snp}; # if the snp is in the ld filtered set already ignore it
            push @snps_filtered, $snp; # if thisis the first time it is seen, add it to the filtered snps, and remove anything in LD with it
            next if $block =~ /NONE/; # nothing is in LD
            my (@block) = split (/\|/, $block);
            foreach my $ldsnp (@block){
                if (exists $snps{$ldsnp}) {
                    $ld_excluded{$ldsnp} = $snp; #Add to the excluded snps, if itis in an LD block with the current snp, and it its one of the test snps.
                    say "$ldsnp excluded for LD at >= $r2 with $snp";
                }
            }
        }
    }

    return (\%ld_excluded, @snps_filtered);#note that if a SNP doesn't exist in the ld file it will be rejected regardless, may need to add these back
}


=head2 get_cells

Read the correct cell list based on data (erc -encode). Also gets the tissue names for the cells.

=cut

sub get_cells{
    # read the correct cell list based on data (erc -encode). Also gets the tissue names for the cells.
    my ($data, $dbh) = @_;
    my $table = "cells_".$data;
    # Check that the table exists in the DB (note, some magic here that might be SQLite-specific)
    my @tables = grep {/^cells_/} map {$_ =~ s/"//g; $_ =~ s/^main\.//; $_} $dbh->tables();
    if (!grep {/$table/} @tables) {
        die "The database does not contain information for the data background provided.\n";
    }
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

=head2 assign

Assign any maf, gc, tss values to the percentile bins

=cut

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

=head1 AUTHOR

Ian Dunham, C<< <dunham at ebi.ac.uk> >>

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Forge


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2014 Ian Dunham.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; version 2 dated June, 1991 or at your option
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available in the source tree;
if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA


=cut

1;
