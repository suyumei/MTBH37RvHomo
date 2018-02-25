#!/usr/bin/perl
#use strict;
use IO::File;
use IO::Handle;
my $goadb="..\\goa\\";
SplitGOA();
print "finished\n";

sub SplitGOA
{
    my %filehandles;
    my $gofile=$goadb."GOA.new.uniprot.proteins";
    open IN, "<$gofile";
    my $fh=IO::Handle->new_from_fd(fileno IN,r);
    while(my $line=$fh->getline)
    {
       chomp($line);
       my @dataitems=split(/\t/,$line);
       my $proteinid=$dataitems[0];
       print "---$proteinid---\n";
       my $initial=substr($proteinid,0,1);
       if(not exists($filehandles{$initial}))
       {
           my $goafile=$goadb."$initial.goa.part";
           $filehandles{$initial}=new IO::File(">$goafile");
       }
       $filehandles{$initial}->print($line."\n")
    }
    $fh->close();
    close IN;
    foreach my $initial(keys%filehandles)
    {
        $filehandles{$initial}->close();
    }
}

