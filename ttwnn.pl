#!/usr/bin/env perl

=head1 NAME

ttwnn.pl - the tool with no name.

=head1 SYNOPSIS

ttwnn.pl options (-f file) (-snps snplist)

=head1 DESCRIPTION

Analyse a set of SNPs for their overlap with DNase 1 hotspots compared to matched background SNPs. Identified Enrichemnt in DHS by tissue.


=head1 OPTIONS

=over

=item B<data>

Data set to analyse. Either ENCODE data ('encode') or Roadmap Epigenome data ('erc'). erc by default.

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

Defaults to 'gwas'. In both cases SNPs have to be on the arrays AND in the 1000 genomes phase 1 integrated call data set at phase1/analysis_results/integrated_call_sets

=item B<label>

Supply a label that you want to use for the plotting titles, and filenames

=item B<f>

Supply the name of a file containing a list of SNPs currently in format chr\tbeg\tend\trsid\tpval. If not supplied the analysis is performed either on snps provided as rsids in a comma separated list through the snps option or on a set of data from a gwas study on Pulmonary_function (http://www.ncbi.nlm.nih.gov/pubmed/21946350, http://www.ncbi.nlm.nih.gov/pubmed/20010835 and http://www.ncbi.nlm.nih.gov/pubmed/20010834)

=item B<snps>

Can provide the snps as a comma separated list.

=item B<format>

bed  = File given is a bed file of locations (chr\tbeg\tend) aka Personal Genome SNP format.  bed format should be 0 based and the chromosome should be given as chrN. Hoever will also accept chomosomes as just N (ensembl) and 1-based format where beg and end are the same

vcf = File given is a vcf file

tabix = File contains SNPs in tabix format.

ian = 1-based chr\tbeg\tend\trsid\tpval

=back

=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http://www.ensembl.org/info/about/code_licence.html for details.

=head1 AUTHOR

Ian Dunham, EBI

=head1 CONTACT

Ian Dunham <dunham@ebi.ac.uk>

=cut

use 5.010;
use warnings;
use DBI;
use Sort::Naturally;
use Cwd;
use Storable;
use Getopt::Long;
use File::Basename;
use Config::IniFiles;


my $cwd = getcwd;

my ($bkgd, $data, $label, $file, $format, @snplist);

GetOptions (
    'data=s'    => \$data,
    'bkgd=s'      => \$bkgd,
    'label=s'   => \$label,
    'f=s'       => \$file,
    'format=s'  => \$format,
    'snp=s'     => \@snplist,
);

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

my $cfg = Config::IniFiles->new( -file => "$dirname/ttwnn.ini" );
my $datadir = $cfg->val('Files', 'datadir');

# percentile bins for the bkgrd calculations. This is hard coded so there are enough SNPs to choose from, but could later be altered.
my $per = 10;
# number of sets to analyse for bkgrd. Again currently hardcoded to 100
my $reps = 100;

unless (defined $bkgd){
    $bkgd = "gwas";
}

my $dsn = "dbi:SQLite:dbname=" . $datadir . "ttwnn.db";
my $dbh = DBI->connect($dsn, "", "") or die $DBI::errstr;
# snps need to come either from a file or a list
my @snps;

if (defined $file){
    # would be better if this was module or sub :-(
    my $sth = $dbh->prepare("SELECT rsid FROM bits WHERE location = ?");
    open my $fh, "<", $file or die "cannot open file $file : $!";
    if ($format =~ /rsid/){
        while (<$fh>){
            chomp;
            push @snps, $_;
        }
    }
    elsif ($format =~ /ian/){
        while (<$fh>){
            my ($chr, $beg, $end, $rsid, undef) = split "\t", $_;
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
    @snps = qw(rs2865531 rs2395730 rs12914385 rs11168048 rs1529672 rs357394 rs13147758 rs3769124 rs2647044 rs12504628 rs1541374 rs2869967 rs1928168 rs3094548 rs3867498 rs9978142 rs4762767 rs6903823 rs11172113 rs9310995 rs2571445 rs2070600 rs11727189 rs3734729 rs2906966 rs1036429 rs16909898 rs3995090 rs12477314 rs2544527 rs2284746 rs993925 rs2277027 rs1344555 rs1455782 rs2855812 rs2838815 rs11001819 rs12716852 rs2798641 rs4129267 rs7068966 rs12899618 rs153916 rs1551943 rs730532 rs1980057 rs3820928 rs2036527 rs10516526 rs2857595 rs3817928 rs310558 rs808225 rs12447804);
}

if (scalar @snps < 20){
    die "Fewer than 20 SNPs provided. Analysis not run";
}
#Connect to the sqlite database file which contains the tables for each data

my $sth = $dbh->prepare("SELECT * FROM bits WHERE rsid IN (?)");

# get the cell list array and the hash that connects the cells and tissues

my ($cells, $tissues) = get_cells($data, $dbh);

# get the bit strings for the test snps from the database file
my $rows = get_bits(\@snps, $sth);

# unpack the bitstrings and store the overlaps by cell.
my $test = process_bits($rows, $cells, $data);

my @missing;
foreach my $rsid (@snps){
    unless (exists $$test{'SNPS'}{$rsid}){
        push @missing, $rsid;
    }
}
if (scalar @missing > 0) {
    print "The following SNPs have not been analysed\n";
    print join("\n", @missing) . "\n";
}

# only pick background snps matching snps that were found originally.
my @foundsnps = keys %{$$test{'SNPS'}};
print "Test SNPs analysed " . scalar @foundsnps . "\n";

# identify the gc, maf and tss, and then make bkgrd picks
my $picks = match(\%$test, $sth, $bkgd);

# for bgrd set need to get distribution of counts instead
# make a hash of data -> cell -> bkgrd-Set -> overlap counts

my %bkgrd; #this hash is going to store the bkgrd overlaps

foreach my $bkgrd (keys %{$picks}){
    $rows = get_bits(\@{$$picks{$bkgrd}}, $sth);
    unless (scalar @$rows == scalar @foundsnps){
        print "Background " . $bkgrd . " only " . scalar @$rows . " SNPs out of " . scalar @foundsnps . "\n";
    }
    #print Dumper $rows;
    my $result = process_bits($rows, $cells, $data);
    #print Dumper $result;
    foreach my $cell (keys %{$$result{'CELLS'}}){
        push @{$bkgrd{$cell}}, $$result{'CELLS'}{$cell}{'COUNT'}; # accumulate the overlap counts by cell
    }
}
$dbh->disconnect();

#Having got the test overlaps and the bkgd overlaps now calculate Zscores and output the table to be read into R for plotting.
my $time = time();
my $resultsdir = "$cwd/$lab.$time";
mkdir $resultsdir;
my $filename = "$lab.$time.chart.tsv";
open my $ofh, ">", "$resultsdir/$filename" or die "Cannot open $resultsdir/$filename: $!"; #should grab a process number for unique name here
print $ofh join("\t", "Zscore", "Cell", "Tissue", "File", "SNPs", "Number") ."\n";
my $n =1;

foreach my $cell (sort {ncmp($$tissues{$a}{'tissue'},$$tissues{$b}{'tissue'}) || ncmp($a,$b)} @$cells){ # sort by the tissues alphabetically (from $tissues hash values)
    # ultimately want a data frame of names(results)<-c("Zscore", "Cell", "Tissue", "File", "SNPs")
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
    my $snp_string = "";
    $snp_string = join(",", @{$$test{'CELLS'}{$cell}{'SNPS'}}) if defined $$test{'CELLS'}{$cell}{'SNPS'}; # This gives the list of overlapping SNPs for use in the tooltips. If there are a lot of them this can be a little useless
    my ($shortcell, undef) = split('\|', $cell); # undo the concatenation from earlier to deal with identical cell names.
    print $ofh join("\t", $zscore, $shortcell, $$tissues{$cell}{'tissue'}, $$tissues{$cell}{'file'}, $snp_string, $n) . "\n";
    $n++;
}

Chart($filename, $lab, $time, $resultsdir); # basic pdf plot
rChart($filename, $lab, $time, $resultsdir); # rCharts polychart plot
dChart($filename, $lab, $time, $resultsdir); # rCharts Dimple chart
table($filename, $lab, $time, $resultsdir);
#hChart("$lab.chart.tsv", $label);


### Subroutines ###

sub match{
    # identifies the bins that each of the snps lies in, and then picks 100 matching SNPs.
    my $snps = shift; # ahash that contians the found snps
    my $sth = shift;
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
        #$rs is the test snp, $rsid os the matched snp.
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
            push @{$picks{$n}}, $rsid; # each $n array is a set of snps matching the test set/
        }
    }
    return \%picks;
}

sub process_bits{
    my ($rows, $cells, $data) = @_;
    my %test;
    foreach my $row (@{$rows}){

        my ($location, $rsid, $sum, $bit, $maf, $tss, $gc);
        if ($data =~ /erc/){
            ($location, $rsid, undef, undef, $bit, $sum, $maf, $tss, $gc) = split("\t", join("\t", @$row));
        }
        else{
            ($location, $rsid, $bit, $sum, undef, undef, $maf, $tss, $gc) = split("\t", join("\t", @$row));
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
    my ($snps, $sth) = @_; #get an array of snps and the db to look at
    my @results;
    foreach my $snp (@$snps){
        $sth->execute($snp);
        my $result = $sth->fetchall_arrayref();
        foreach my $row (@{$result}){
            push @results, $row;
        }
    }
    $sth->finish();
    return \@results;# return the bitstring line from the database
}

sub fetch_rsid{
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
    my $data = shift;
    my $dbh = shift;
    my $table = join('_', "cells", $data);
    my $sth = $dbh->prepare("SELECT shortcell,tissue,file FROM $table");
    $sth->execute();
    my $ver = $sth->fetchall_arrayref();
    $sth->finish();
    my ($cells, $tissues);
    foreach my $row (@$ver){
        my $cell = shift @$row;
        my $tissue = shift @$row;
        my $file = shift @$row;

        $cell = "$cell|$file"; # Sometimes the same cell is used twice, with a differnt file so need to record separately (e.g. WI-38).
        push @$cells, $cell;
        $$tissues{$cell}{'tissue'} = $tissue; # this is the hash that is used to connect cells and tissues and ultimately provide the sorting order
        $$tissues{$cell}{'file'} = $file;
    }
    #print Dumper %$tissues;
    return ($cells, $tissues); # return
}

sub assign{
    #sub routine to assign any maf, gc, tss values to the percentile bin
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
    #
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
    #
    my $ev = mean(@_);
    my $sum = 0;
    foreach (@_) { $sum += ($_ - $ev)**2 };

    return $sum/($#_+1);
}

# calulates the standard deviation of an array: this is just the sqrt of the var
sub std { sqrt(var(@_)) }

sub Chart{
    # This is the original code using standard R plot to generate a static pdf.
    print "Making static chart.\n";
    my ($filename, $lab, $time, $resultsdir) = @_;
    my $Rdir = $resultsdir;
    my $chart = "$lab.$time.chart.pdf";
    my $rfile = "$Rdir/$lab.$time.chart.R";

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
results\$Class<-cut(results\$Zscore, breaks =c(min(results\$Zscore), 2.58, 3.39, max(results\$Zscore)), labels=FALSE, include.lowest=TRUE) # 99.9 and 99% CIs 1, 2, 3
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
abline(h=-2.58, col=\"lightpink1\") # Z score of 2.58 = 99 % probability
abline(h=2.58, col=\"lightpink1\")
abline(h=3.39, col=\"lightpink1\", lty=2)
abline(h=3.39, col=\"lightpink1\", lty=2)
text(c(-2),2.58+0.14,c(\"99%\"),col=\"lightpink1\",adj=1,cex=0.8)
text(c(-2),-2.58+0.14,c(\"1%\"),col=\"lightpink1\",adj=1,cex=0.8)
text(c(-1),3.39+0.14,c(\"99.9%\"),col=\"lightpink1\",adj=1,cex=0.8)
text(c(-1),-3.39+0.14,c(\"0.1%\"),col=\"lightpink1\",adj=1,cex=0.8)
palette(\"default\")\n";

    foreach my $pos (@lines){
        print $rfh "lines(c($pos,$pos),c(-5,22),lty=6,col=\"#00000070\")\n" unless $pos == 0.5;
    }
    $index = 0;
    foreach my $tissue (@labels){
        print $rfh "text(c(" . $label_pos[$index] . "),ymax,c(\"" . $tissue . "\"),col=\"burlywood3\",adj=1,srt=90,cex=0.8)\n";
        $index++;
    }

    print $rfh "dev.off()\n";

    system "R --no-save --quiet --slave < $rfile";
}

sub rChart{
    print "Making rChart.\n";
    my ($filename, $lab, $time, $resultsdir) = @_;
    my $chart = "$lab.$time.rchart.html";
    my $Rdir = $resultsdir;
    my $rfile = "$Rdir/$lab.$time.rChart.R";
    open my $rcfh, ">", "$rfile";
    print $rcfh "setwd(\"$Rdir\")
results<-read.table(\"$filename\", header = TRUE, sep=\"\t\")
results\$Colour<- 0 + (results\$Zscore < 3.39) + (results\$Zscore < 2.58)  # 99.9 and 99% CIs
require(rCharts)
r1 <- rPlot(Zscore ~ Cell, data=results, color=\"bin(Colour, 0.25)\", type=\"point\", tooltip = \"function(item){ return (item.Zscore + '\\n' + item.Cell + '\\n' + item.Tissue + '\\n' + item.File + '\\n' + item.SNPs + '\\n' )}\")
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
    print "Making dChart.\n";
    my ($filename, $lab, $time, $resultsdir) = @_;
    my $chart = "$lab.$time.dchart.html";
    my $Rdir = $resultsdir;
    my $rfile = "$Rdir/$lab.$time.dChart.R";
   open my $rcfh, ">", "$rfile";
    print $rcfh "setwd(\"$Rdir\")
results<-read.table(\"$filename\", header = TRUE, sep=\"\t\")
results\$Class<-cut(results\$Zscore, breaks =c(min(results\$Zscore), 2.58, 3.39, max(results\$Zscore)), labels=FALSE, include.lowest=TRUE) # 99.9 and 99% CIs 1, 2, 3
require(rCharts)
d1 <- dPlot(
  y = \"Zscore\",
  x = c(\"Cell\", \"Tissue\", \"SNPs\", \"Number\"),
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
    print "Making Table.\n";
    my ($filename, $lab, $time, $resultsdir) = @_;
    my $chart = "$lab.$time.table.html";
    my $Rdir = $resultsdir;
    my $rfile = "$Rdir/$lab.$time.table.R";
    open my $rcfh, ">", "$rfile";
    print $rcfh "setwd(\"$Rdir\")
    data<-read.table(\"$filename\", header = TRUE, sep=\"\t\")
    results<-data.frame(data\$Cell, data\$Tissue, data\$Zscore, data\$SNPs)
    names(results)<-c(\"Cell\", \"Tissue\", \"Zscore\", \"SNPs\")
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
#results\$Class<-cut(results\$Zscore, breaks =c(min(results\$Zscore), 2.58, 3.39, max(results\$Zscore)), labels=FALSE, include.lowest=TRUE) # 99.9 and 99% CIs 1, 2, 3
#require(rCharts)
#h1 <- hPlot(Zscore ~ Number, data=results, type=\"scatter\", radius=5)
#h1\$addParams(width = 2000, height=600, title=list(\"$label overlaps with $data DHS\"))
#h1\$xAxis(title = list(text = \"Cell\"), labels = list(rotation=-90, align=\"right\"), categories = results\$Cell\)
#h1\$save('$chart', cdn = F)";
#
#system "R --no-save --quiet --slave < $rfile";
#}
