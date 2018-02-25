#!/usr/bin/perl
use IO::File;
use IO::Handle;

my $fastafile=$ARGV[0];
#my $fastafile='test';
my $datadir="..\\data\\";
my $featuredir="..\\featurematrix\\";
my $resultdir="..\\results\\";
my $topic='MTBH37RvHomo';

GenerateResults();
print "finished\n";


sub GenerateResults
{
    my $no=0;
    my %predictionset;
    my $predictionsetfile="$datadir$fastafile";
    open IN, "<$predictionsetfile";
    my $fh=IO::Handle->new_from_fd(fileno IN,r);
    while(my $line=$fh->getline)
    {
       chomp($line);
       $no++;
       my @tt=split(/\t/,$line);
       $predictionset{$no}=$tt[0]."\t".$tt[1];
    }
    close IN;
    $fh->close();

    my $resultfile=$resultdir."$fastafile.results.txt";
    open my $resultfilehandle, ">$resultfile" ;

    my $no=0;
    print "loading prediction result file:...\n";
    my $predictionresultfile=$featuredir."$fastafile.predict.final.decvalues.txt";
    open IN, "<$predictionresultfile";
    my $fh=IO::Handle->new_from_fd(fileno IN,r);
    while(my $line=$fh->getline)
    {
       chomp($line);
       $no++;
       if($line>0)
       {
          print $resultfilehandle  $predictionset{$no}."\t".$line."\tinteraction\n";

       }
       elsif($line<0)
       {
          print $resultfilehandle  $predictionset{$no}."\t".$line."\tnon-interaction\n";
       }
       else
       {
          print $resultfilehandle  $predictionset{$no}."\t".$line."\tundetermined\n";
       }
    }
    close IN;
    $fh->close();
    close $resultfilehandle;
}