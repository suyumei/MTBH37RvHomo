First, unzip all the zipped files and place the unzepped files directly under the subdirectory 'goa'.
  -blast/data:swissprot.zip,swissprot.psq.zip
  -goa:A.goa.zip,D.goa.zip,F.goa.zip,G.goa.zip,Q.goa.zip

WARNING: the source codes would fail to work if the requested files (blast/data:swissprot,swissprot.psq;goa:A.goa.part,D.goa.part,F.goa.part,G.goa.part,Q.goa.part) are not found in the subdirectories 'blast/data' and 'goa'.
These files have to be zipped due to github 100MB limit of files.

1. Workflows:
  -----------------------
(1) Data preparation: prepare the query file and its corresponding fasta file, and place the files into the subdirectory 'data'. 
-Each protein should be given as the format genesymble|Uniprot accession (e.g. fadD5|O07411) (see test);
-Each entry in the fasta file (query file name suffixed with .fasta) should be given as the format (see test.fasta)
  >uniprot accession 
  sequence

(2) Run PSIBlastlargedatadetail.pl in the subdirectory 'code' to obtain homologs.

(3) Run GenerateFeatureMatrix.pl in the subdirectory 'code' to obtain feature representation. 
-This step is time-consuming, advanced users are encouraged to deploy local GOA database.
-Step (2) and (3) could also be called in the script Main.m, but no progress information would appear in Matlab.

(4) Run the script Main.m in the subdirectory 'code' to get results, the results can be found in the subdirectory 'results' (e.g. test.results.txt)

2. Subdirectories
  -----------------------
-data: query file and fasta file;
-code: source code (matlab and perl script);
-featurematrix: GO term & index, model file, feature vectors. The file initialed with MTBH37RvHomo should not be deleted;
-goa: goa files. To accelerate goa searching, the go annotations are arranged in the order of the initial character of uniprot accession (A-Z), advanced users are encouraged to deploy local GOA database;
-results: prediction results;
-blast: psiblast tool

3. Supported platforms 
  -----------------------
MS Windows 32-bit/64-bit, matlab 2015, activeperl 5.24.1

4. Demo
  -----------------------
Run codes as specified in the workflows given the examples in the file test.

5. Applicability
  -----------------------
Directly replace the protein pairs of the file test and the corresponding fasta information of the file test.fasta in the subdirectory 'data'.