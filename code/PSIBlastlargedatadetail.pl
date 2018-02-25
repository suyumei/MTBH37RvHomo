#!/usr/bin/perl
#use strict;

use IO::File;
use IO::Handle;

my $fastafile=$ARGV[0];
my $fastafile='test';
my $topic='MTBH37RvHomo';
my $datadir="..\\data\\";
my $blastbindir="..\\blast\\bin\\";
my $blastdatadir="..\\blast\\data\\";
my $featuredir="..\\featurematrix\\";
my $tmp_file=$featuredir."tmp";
my @fastafiles;
push(@fastafiles,$fastafile);
my %blastaccessions;
my %downloadedproteins;
print "psi-blasting...\n";
PsiblastHomology();
print "finished...\n";

#read fasta sequence data
sub PsiblastHomology
{
    print "Reading total fasta data...\n";

    for(my $i=0;$i<@fastafiles;$i++)
    {
          my $fasta_file=$fastafiles[$i];
          my $data_file=$datadir."$fasta_file.fasta";
          print $data_file."\n";

          $downloadedfile=$featuredir."$fasta_file.downloaded";
          open IN, "<$downloadedfile";
           my @querylines = <IN>;
           close IN;

           for(my $i=0;$i<@querylines;$i++)
           {
               $lastaccession=$querylines[$i];
               chomp($lastaccession);
               $downloadedproteins{$lastaccession}="";
           }

          open IN, "<$data_file" ;
          my @querylines = <IN>;
          close IN;
          my $queryname="";
          my $queryseq="";

          for(my $j=0;$j<@querylines;$j++)
          {
              my $query=$querylines[$j];
              chomp($query);
              if($query=~m/^>/)
              {

                 if($queryseq ne "")
                 {
                     PSIBLAST($fasta_file,$queryname,$queryseq);
                 }
                 $queryname=$query;
                 $queryname=~s/>//;
                 $queryseq="";
              }
              else
              {
                  if($queryseq eq "")
                  {
                      $queryseq= $query;
                  }
                  else
                  {
                      $queryseq= $queryseq.$query;
                  }
              }
          }
          #last query
          if($queryname ne "")#complete read last sequence
          {
              PSIBLAST($fasta_file,$queryname,$queryseq);
          }

    }
}
sub PSIBLAST
{
    my $fasta_file=shift;
    my $queryname=shift;
    my $queryseq=shift;
    if(not exists($downloadedproteins{$queryname}))
    {
        print "Psi-blasting ".$queryname."...\n";
        open $tmp_filehandle, ">$tmp_file";
        print $tmp_filehandle ">".$queryname."\n".$queryseq;
        close $tmp_filehandle;
        my @homology;
        GetHomologs($queryname,$tmp_file,\@homology);

        my $homologyaccessions=$queryname;
        for($k=0;$k<@homology;$k++)
        {
            $homologyaccessions=$homologyaccessions.";".$homology[$k];
        }
        my $homologyfile=$featuredir."$fasta_file.homology";
        open $homologyfilehandle,">>$homologyfile";
        print $homologyfilehandle $homologyaccessions."\n";
        close $homologyfilehandle;

        my $downloadedfile=$featuredir."$fasta_file.downloaded";
        open my $downloadedfilehandle,">>$downloadedfile";
        print $downloadedfilehandle $queryname."\n";
        close $downloadedfilehandle;
        $downloadedproteins{$queryname}="";

    }
}
sub GetHomologs
{
    my $queryname;
    my $tmp_file;
    my $homology;
    ($queryname,$tmp_file,$homology)=@_;
    $algnfile=$tmp_file.".algn";
    if($^O eq "MSWin32")
    {
         $cmdline=sprintf("%sblastall.exe -p blastp  -m 8 -d  %sswissprot -i $tmp_file -o $algnfile",$blastbindir,$blastdatadir);
    }
    else
    {
         $cmdline="blastp -d uniprot -i $tmp_file -o $algnfile -m 6 -j 1 -e 0.0001 -h 0.002";
    }
    system($cmdline);
    if (not -e $algnfile)
    {
        print "failed to blast out algn file\n";
    }

    open my $IN, "<$algnfile" or die "Cannot open input file $algnfile for reading\n";
    my @lines = <$IN>;
    close $IN;
    my $counter=0;
    my %homologs;
    for(my $i=0;$i<@lines;$i++)
    {
        $temp=$lines[$i];
        chomp($temp);
        @t=split(/\t/,$temp);
        my $name=$t[1];
        my $identity=$t[2];
        my $evalue=$t[10];
        my @tt=split(/\|/,$name);
        my $homologyname=$tt[3];
        $homologyname=(split(/\./,$homologyname))[0];
        my $species=$tt[4];

        if(not exists($homologs{$homologyname}))
        {
            $homologs{$homologyname}="";
            $homology->[$counter]=$homologyname.":".$identity.":".$evalue.":".$species;
            $counter++;
        }
    }
    if($counter==0)
    {
        print "no homology found!\n"
    }
}