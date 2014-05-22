package Forge::Plot;

use 5.010;
use strict;
use warnings FATAL => 'all';
use Sort::Naturally;

=head1 NAME

Forge::Plot - Plotting utilities for Forge

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

our (@ISA, @EXPORT, @EXPORT_OK);
use Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(Chart dChart table); # Symbols to export by default
@EXPORT_OK = qw(rChart);

=head1 SYNOPSIS

Provise plotting utilities for differnt plots to Forge

=head1 EXPORT

Chart
dChart
table

=head1 SUBROUTINES/METHODS

=head2 Chart

This is the original code using standard R plot to generate a static pdf.

=cut


sub Chart{
    print "Making static chart.\n";
    my ($filename, $lab, $resultsdir, $tissues, $cells, $label, $t1, $t2, $data) = @_;
    my $Rdir = $resultsdir;
    my $chart = "$lab.chart.pdf";
    my $rfile = "$Rdir/$lab.chart.R";
    #set some colors
    my ($sig, $msig, $ns, $abline, $tline) = qw(red palevioletred1 steelblue3 lightpink1 burlywood3); #alternate msig = pink2

    # make plot, first calculate where dividing lines are:
    my (@lines, @label_pos, @labels, @tissue_txt);
    my $n = 1;
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
    my $index = 0;
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
#results\$Class<-cut(results\$Pvalue, breaks =c(min(results\$Pvalue), $t1, $t2, max(results\$Pvalue)), labels=FALSE, include.lowest=TRUE) # 99.9 and 99% CIs 1, 2, 3
$t1 = sprintf("%.2f", $t1);
$t2 = sprintf("%.2f", $t2);
    print $rfh "setwd(\"$Rdir\")
results<-read.table(\"$filename\",header=TRUE,sep=\"\t\")
results\$Class<-cut(results\$Pvalue, breaks =c(min(results\$Pvalue), $t1, $t2, max(results\$Pvalue)), labels=FALSE, include.lowest=TRUE)
pdf(\"$chart\", width=22.4, height=7)
palette(c(\"$ns\",\"$msig\",\"$sig\"))
#ymin1 = min(results\$Pvalue, na.rm=TRUE)*1.1
ymin1= min(results\$Pvalue, na.rm=TRUE)*1.1
ymax1 = max(results\$Pvalue, na.rm=TRUE)*1.1
ymax = max(c(abs(ymin1),ymax1))
ymin = -ymin1
par(mar=c(9,4,3,1)+0.1)
plot(results\$Pvalue,ylab=\"-log10 binomial P\",xlab=\"\",main=\"SNPs in DNase1 sites (probably TF sites) in cell lines for $data $label\",ylim=c(ymin,ymax), las=2, las=2, pch=19,col=results\$Class, xaxt='n')
axis(1, seq(1,length(results\$Cell)),labels=results\$Cell, las=2, cex.axis=0.7)
mtext(1,text=\"Cell\",line=7,cex=1.2)
abline(h=$t1, col=\"$abline\")
#abline(h=-$t2, col=\"$abline\", lty=2)
abline(h=$t2, col=\"$abline\", lty=2)
text(c(-0.5),$t1+0.2,c(\"P = $t1\"),col=\"$abline\",adj=1,cex=0.8)
#text(c(-0.5),-$t1+0.16,c(\"1%\"),col=\"$abline\",adj=1,cex=0.8)
text(c(-0.5),$t2+0.2,c(\"P = $t2\"),col=\"$abline\",adj=1,cex=0.8)
#text(c(-0.5),-$t2+0.16,c(\"0.1%\"),col=\"$abline\",adj=1,cex=0.8)
palette(\"default\")\n";

    foreach my $pos (@lines){
        print $rfh "lines(c($pos,$pos),c(-22,22),lty=6,col=\"#00000070\")\n" unless $pos == 0.5;
    }
    $index = 0;
    foreach my $tissue (@labels){
        print $rfh "text(c(" . $label_pos[$index] . "),ymax,c(\"" . $tissue . "\"),col=\"$tline\",adj=1,srt=90,cex=0.8)\n";
        $index++;
    }
    print $rfh "dev.off()\n";
#run the R code
    system "R --no-save --quiet --slave < $rfile";
}

=head1 SUBROUTINES/METHODS

=head2 rChart

Makes a polycharts plot : note X axis labelling is problematical

=cut

sub rChart{

    print "Making rChart.\n";
    my ($filename, $lab, $resultsdir, $data, $label, $t1, $t2) = @_;
    my $chart = "$lab.rchart.html";
    my $Rdir = $resultsdir;
    my $rfile = "$Rdir/$lab.rChart.R";
    open my $rcfh, ">", "$rfile";
    print $rcfh "setwd(\"$Rdir\")
results<-read.table(\"$filename\", header = TRUE, sep=\"\t\")
results\$Colour<- 0 + (results\$Pvalue < $t2) + (results\$Pvalue < $t1)  # 99.9 and 99% CIs
require(rCharts)
r1 <- rPlot(Pvalue ~ Cell, data=results, color=\"bin(Colour, 0.25)\", type=\"point\", tooltip = \"function(item){ return (item.Pvalue + '\\\\n' + item.Cell + '\\\\n' + item.Tissue + '\\\\n' + item.File + '\\\\n' + item.SNPs + '\\\\n' + item.Accession + '\\\\n')}\")
#r1\$guides(color=list(scale = list(type = \'gradient\', lower = \'\#CCC\', upper = \'\#000\'))) # optional code to make a grey scale
r1\$addParams(width = 2000, height=600, title=\"$label overlaps with $data DHS\")
ymin1 = min(results\$Pvalue, na.rm=TRUE)*1.2
ymax1 = max(results\$Pvalue, na.rm=TRUE)*1.2
ymax = max(c(abs(ymin1),ymax1))
ymin = -ymax
r1\$guides(x = list(numticks = length(unique(results\$Cell)), levels=results\$Cell), y = list(min = ymin, max = ymax))
r1\$save('$chart', cdn = F)
##r1\$show() #makes a temp file\n";

system "R --no-save --quiet --slave < $rfile";
}

=head1 SUBROUTINES/METHODS

=head2 dChart

Make dimple interactive chart.

=cut

sub dChart{

    print "Making dChart.\n";
    my ($filename, $lab, $resultsdir, $data, $label, $t1, $t2) = @_;
    my $chart = "$lab.dchart.html";
    my $Rdir = $resultsdir;
    my $rfile = "$Rdir/$lab.dChart.R";
   open my $rcfh, ">", "$rfile";
    print $rcfh "setwd(\"$Rdir\")
results<-read.table(\"$filename\", header = TRUE, sep=\"\t\")
results\$Class<-cut(results\$Pvalue, breaks =c(min(results\$Pvalue), $t1, $t2, max(results\$Pvalue)), labels=FALSE, include.lowest=TRUE)
d1 <- dPlot(
  y = \"Pvalue\",
  x = c(\"Cell\", \"Tissue\", \"SNPs\", \"Number\", \"Accession\", \"Pvalue\", \"Zscore\"),
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

=head1 SUBROUTINES/METHODS

=head2 table

=cut

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
    results<-data.frame(data\$Cell, data\$Tissue, data\$Accession, data\$Pvalue, data\$Zscore, data\$SNPs)
    names(results)<-c(\"Cell\", \"Tissue\", \"Accession\", \"Pvalue\", \"Zscore\", \"SNPs\")
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

=head1 SUBROUTINES/METHODS

=head2 hChart

highcharts interface has problems with plotting - not used

=cut

#
#sub hChart{
#    print "Making hChart.\n";
#    my ($filename, $label) = @_;
#    my $chart = $lab . ".hchart.html";
#    my $Rdir = $cwd;
#    my $rfile = $lab . ".hChart.R";
#    open my $rcfh, ">", "$Rdir/$rfile";
#    print $rcfh "setwd(\"$Rdir\")
#results<-read.table(\"$filename\", header = TRUE, sep=\"\t\")
#results\$Class<-cut(results\$Pvalue, breaks =c(min(results\$Pvalue), $t1, $t2, max(results\$Pvalue)), labels=FALSE, include.lowest=TRUE) # 99.9 and 99% CIs 1, 2, 3
#require(rCharts)
#h1 <- hPlot(Pvalue ~ Number, data=results, type=\"scatter\", radius=5)
#h1\$addParams(width = 2000, height=600, title=list(\"$label overlaps with $data DHS\"))
#h1\$xAxis(title = list(text = \"Cell\"), labels = list(rotation=-90, align=\"right\"), categories = results\$Cell\)
#h1\$save('$chart', cdn = F)";
#
#system "R --no-save --quiet --slave < $rfile";
#}



=head1 AUTHOR

Ian Dunham, C<< <dunham at ebi.ac.uk> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-forge at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Forge>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Forge::Plot


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

1; # End of Plot
