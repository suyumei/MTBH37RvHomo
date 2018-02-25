clear;
params.topic='MTBH37RvHomo';
params.relativePath='..\\featurematrix\\';
params.queryfile='test';
params.vectortypes={'.feature.vector.txt'};

% disp('PSLBLASTing...');
% perl_cmd = sprintf('!perl PSIBlastlargedatadetail.pl %s',params.queryfile);
% eval(perl_cmd); 
% disp ('searching GOA files for feature representation,this step is time-consuming (advanced users are encouraged to deploy local GOA database,or run GenerateFeatureMatrix.pl independently), please wait...');
% strcmd = sprintf('perl GenerateFeatureMatrix.pl %s',params.queryfile);
% eval(strcmd);     

disp ('FormatFeatureMatrix.pl...');
strcmd = sprintf('perl FormatFeatureMatrix.pl %s',params.queryfile);
eval(strcmd);
PredictMain(params);
disp('finished!');
