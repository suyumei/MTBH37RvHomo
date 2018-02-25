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
my $goadb="..\\goa\\";

my @proteinsIDs;
my %initials;
my %totalProteinIDs;
my %proteinhomologies;
my %proteinGOterms;
my %Fgotermindex;
my %Cgotermindex;
my %Pgotermindex;
my $Fgotermcounter=0;
my $Cgotermcounter=0;
my $Pgotermcounter=0;
my @filehandles;
my @fastafiles;
push(@fastafiles,$fastafile);

print "reading homology data...\n";
GetProteinhomologies();
print "searching go terms...\n";
LoadGOTermIndex();
GetQueryProteinGOTerms();
print "generating go term vector...\n";
GenerateQueryProteinGOFeatureMatrix();
print "finished...\n";

sub GetProteinhomologies
{
    #for(my $i=0;$i<@fastafiles;$i++)
    {
        my $homologyfile=$featuredir."$fastafile.homology";
        open IN, "<$homologyfile" ;
        my @querylines = <IN>;
        close IN;
        for(my $j=0;$j<@querylines;$j++)
        {
            my $query=$querylines[$j];
            chomp($query);
            @t=split(/;/,$query);
            my $proteinID=$t[0];
            $proteinhomologies{$proteinID}=$query;
            for(my $k=0;$k<@t;$k++)
            {
                my @tt=split(/:/,$t[$k]);
                my $p=$tt[0];
                my $initial=substr($p,0,1);
                $initials{$initial}='';
                if(not exists($totalProteinIDs{$p}))
                {
                   $totalProteinIDs{$p}="true";
                }
            }
        }
    }
}

sub LoadGOTermIndex
{
    my $gofile=$featuredir.$topic.".target.homolog.go.terms";
    open IN, "<$gofile";
    my $fh=IO::Handle->new_from_fd(fileno IN,r);
    while(my $line=$fh->getline)
    {
          chomp($line);
          my @dataitems=split(/\t/,$line);
          my $proteinid=$dataitems[0];
          my $goterms=$dataitems[1];
          $proteinGOterms{$proteinid}=$goterms;
    }
    $fh->close();
    close IN;

    my $goindexfile=$featuredir.$topic.".go.terms.index.C";
    open IN, "<$goindexfile";
    my $fh=IO::Handle->new_from_fd(fileno IN,r);
    while(my $line=$fh->getline)
    {
          chomp($line);
          my @dataitems=split(/\t/,$line);
          my $goterm=$dataitems[0];
          my $index=$dataitems[1];
          $Cgotermindex{$goterm}=$index;
          $Cgotermcounter++;
    }
    $fh->close();
    close IN;

    my $goindexfile=$featuredir.$topic.".go.terms.index.F";
    open IN, "<$goindexfile";
    my $fh=IO::Handle->new_from_fd(fileno IN,r);
    while(my $line=$fh->getline)
    {
          chomp($line);
          my @dataitems=split(/\t/,$line);
          my $goterm=$dataitems[0];
          my $index=$dataitems[1];
          $Fgotermindex{$goterm}=$index;
          $Fgotermcounter++;
    }
    $fh->close();
    close IN;

    my $goindexfile=$featuredir.$topic.".go.terms.index.P";
    open IN, "<$goindexfile";
    my $fh=IO::Handle->new_from_fd(fileno IN,r);
    while(my $line=$fh->getline)
    {
          chomp($line);
          my @dataitems=split(/\t/,$line);
          my $goterm=$dataitems[0];
          my $index=$dataitems[1];
          $Pgotermindex{$goterm}=$index;
          $Pgotermcounter++;
    }
    $fh->close();
    close IN;
}

sub GetQueryProteinGOTerms
{
    foreach my $c(sort{ord($a)<=>ord($b)}keys%initials)
    {
       my $gofile=$goadb."$c.goa.part";
       open IN, "<$gofile";
       my $fh=IO::Handle->new_from_fd(fileno IN,r);
       while(my $line=$fh->getline)
       {
          chomp($line);
          my @dataitems=split(/\t/,$line);
          my $proteinid=$dataitems[0];
          print "---$proteinid---\n";
          next if(not exists($totalProteinIDs{$proteinid}));
          $proteinGOterms{$proteinid}=$dataitems[1];
          print "---$proteinid found---\n";
       }
       $fh->close();
       close IN;
    }
}

sub GenerateQueryProteinGOFeatureMatrix
{
    for(my $i=0;$i<@fastafiles;$i++)
    {
        my $fasta_file=$fastafiles[$i];
        $targethomologinstancefilename=$featuredir.$fastafiles[$i].".feature.vectors.txt";
        $filehandles[$i]=new IO::File(">$targethomologinstancefilename");
    }
    $numofclass=@fastafiles;
    my @classbagcounter=(0)x$numofclass;
    my @bagno=(0)x$numofclass;
    for(my $i=0;$i<@fastafiles;$i++)
    {
          my $fasta_file=$fastafiles[$i];
          $data_file="$datadir$fasta_file";
          open IN, "<$data_file" ;
          my @querylines = <IN>;
          close IN;
          for(my $j=0;$j<@querylines;$j++)
          {
              my $query=$querylines[$j];
              chomp($query);
              my @t=split(/\s+/,$query);
              my $first=$t[0];
              my $second=$t[1];
              $first=(split(/\|/,$first))[1];#gene|protein
              $second=(split(/\|/,$second))[1];
              print "processing ".($j+1).":{$first,$second}...\n";
              $bagno[$i]++;

              my $targetgofeaturevector;
              my $homologgofeaturevector;

              my $targetVector=$bagno[$i];
              my $homologVector=$bagno[$i];


              getfeaturevector($first,\$targetgofeaturevector,\$homologgofeaturevector);
              my $result;
              my @firsttargetvectors=split(/;/,$targetgofeaturevector);
              my @firsthomologvectors=split(/;/,$homologgofeaturevector);

              getfeaturevector($second,\$targetgofeaturevector,\$homologgofeaturevector);
              my $result;
              my @secondtargetvectors=split(/;/,$targetgofeaturevector);
              my @secondhomologvectors=split(/;/,$homologgofeaturevector);

              my $r1,$r2;
              isAnnotated(\@firsttargetvectors,\$r1);
              isAnnotated(\@secondtargetvectors,\$r2);
              if($r1==1 && $r2==1)
              {
                  $targetVector=$targetVector." 1";
              }
              else
              {
                  $targetVector=$targetVector." 0";
              }

              my $r1,$r2;
              isAnnotated(\@firsthomologvectors,\$r1);
              isAnnotated(\@secondhomologvectors,\$r2);
              if($r1==1 && $r2==1)
              {
                  $homologVector=$homologVector." 1";
              }
              else
              {
                  $homologVector=$homologVector." 0";
              }


              my @result1;
              processfeaturevector($firsttargetvectors[0],$secondtargetvectors[0],\@result1);
              if(@result1==0)
              {
                  @result1=(0)x$Fgotermcounter;
              }

              my $targetFVector=' '.join(' ',@result1);# F


              my @result2;
              processfeaturevector($firsttargetvectors[1],$secondtargetvectors[1],\@result2);
              if(@result2==0)
              {
                  @result2=(0)x$Cgotermcounter;
              }
              my $targetCVector=' '.join(' ',@result2);# C

              my @result3;
              processfeaturevector($firsttargetvectors[2],$secondtargetvectors[2],\@result3);
              if(@result3==0)
              {
                  @result3=(0)x$Pgotermcounter;
              }
              my $targetPVector=' '.join(' ',@result3);# P
              $targetVector=$targetVector.$targetFVector.$targetCVector.$targetPVector;

              my @result4;
              processfeaturevector($firsthomologvectors[0],$secondhomologvectors[0],\@result4);
              if(@result4==0)
              {
                  @result4=(0)x$Fgotermcounter;
              }

              my $homologFVector=' '.join(' ',@result4);# F

              my @result5;
              processfeaturevector($firsthomologvectors[1],$secondhomologvectors[1],\@result5);
              if(@result5==0)
              {
                  @result5=(0)x$Cgotermcounter;
              }
              my $homologCVector=' '.join(' ',@result5);# C

              my @result6;
              processfeaturevector($firsthomologvectors[2],$secondhomologvectors[2],\@result6);
              if(@result6==0)
              {
                  @result6=(0)x$Pgotermcounter;
              }
              my $homologPVector=' '.join(' ',@result6); # P
              $homologVector=$homologVector.$homologFVector.$homologCVector.$homologPVector;

              $filehandles[$i]->print($targetVector."\n");
              $filehandles[$i]->print($homologVector."\n");
          }
          $filehandles[$i]->close();
    }

}


sub isAnnotated
{
    my $t;
    my $r;
    ($t,$r)=@_;
    my @t1=split(/ /,$t->[0]);
    my @t2=split(/ /,$t->[1]);
    my @t3=split(/ /,$t->[2]);
    if((grep {$_>0} @t1)||(grep {$_>0} @t2)||(grep {$_>0} @t3))
    {
       $$r=1;#annotated protein:the 2nd col is 1
    }
    else
    {
       $$r=0;#completely unannotated protein:the 2nd col is 0
    }

}

sub processfeaturevector
{
    my $first;
    my $second;
    my $result;
    ($first,$second,$result)=@_;
    my @firstpart=split(/ /,$first);
    my @secondpart=split(/ /,$second);
    for(my $i=0;$i<@firstpart;$i++)
    {
        my $firstvalue=$firstpart[$i];
        my $secondvalue=$secondpart[$i];
        if(($firstvalue==1) && ($secondvalue==1))
        {
            $result->[$i]=2;
        }
        elsif(($firstvalue==1) && ($secondvalue==0))
        {
            $result->[$i]=1;
        }
        elsif(($firstvalue==0) && ($secondvalue==1))
        {
            $result->[$i]=1;
        }
        else
        {
            $result->[$i]=0;
        }
    }
}
sub getfeaturevector
{
    my $proteinid;
    my $targetgofeaturevector;
    my $homologgofeaturevector;
    ($proteinid,$targetgofeaturevector,$homologgofeaturevector)=@_;
    my $homologies=$proteinhomologies{$proteinid};
    my @tmp=split(/;/,$homologies);

    my @Fgotermfeaturevector=(0)x$Fgotermcounter;
    my @Cgotermfeaturevector=(0)x$Cgotermcounter;
    my @Pgotermfeaturevector=(0)x$Pgotermcounter;

    for(my $j=0;$j<@tmp;$j++)
    {
        my @homologinfo=split(/:/,$tmp[$j]);
        my $homology=$homologinfo[0];
        my $seqsim=$homologinfo[1];
        my $significance=$homologinfo[2];
        next if(($j>0) &&($homology eq $proteinid));

        my $instacneno=$j;
        my $goterms=$proteinGOterms{$homology};
        my @tt=split(/;/,$goterms);
        for(my $k=0;$k<@tt;$k++)
        {
            my @w=split(/\+/,$tt[$k]);
            my $goterm=$w[0];
            my $type=$w[1];
            SWITCH:
            {
                   $type=~m/F/ and do
                   {
                          if(exists($Fgotermindex{$goterm}))
                          {
                              my $index=$Fgotermindex{$goterm};
                              $Fgotermfeaturevector[$index]=1;
                          }
                          last;
                   };
                   $type=~m/C/ and do
                   {
                          if(exists($Cgotermindex{$goterm}))
                          {
                              my $index=$Cgotermindex{$goterm};
                              $Cgotermfeaturevector[$index]=1;
                          }
                          last;
                   };
                   $type=~m/P/ and do
                   {
                          if(exists($Pgotermindex{$goterm}))
                          {
                               my $index=$Pgotermindex{$goterm};
                               $Pgotermfeaturevector[$index]=1;
                          }
                          last;
                   };

            }
        }
        if($j==0)
        {
                $$targetgofeaturevector=join(' ',@Fgotermfeaturevector).';'.join(' ',@Cgotermfeaturevector).';'.join(' ',@Pgotermfeaturevector);
                @Fgotermfeaturevector=(0)x$Fgotermcounter;
                @Cgotermfeaturevector=(0)x$Cgotermcounter;
                @Pgotermfeaturevector=(0)x$Pgotermcounter;

        }
    }
    $$homologgofeaturevector=join(' ',@Fgotermfeaturevector).';'.join(' ',@Cgotermfeaturevector).';'.join(' ',@Pgotermfeaturevector);
}