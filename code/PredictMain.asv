function PredictMain(params)    
    params.Cs=2.^[17];     
    for i=1:size(params.Cs,2)
        params.C=params.Cs(i);
        disp ('predicting...');
        Predict(params);        
    end
    disp ('GetPredictResults.pl...');
    strcmd = sprintf('perl GetPredictResults.pl %s',params.queryfile);
    eval(strcmd);
    disp ('finished!');
    
    function Predict(params)           
        modelfile=sprintf('%s%s.train.multi.labels.model',params.relativePath,params.topic);        
        predictdatafile=sprintf('%s%s.feature.vectors.txt.predict.multi.labels',params.relativePath,params.queryfile); 
        disp ('predicting...'); 
        cmdstr=sprintf('predict.exe -b 1 %s %s %s.predict_labels',predictdatafile,modelfile,predictdatafile);             
        system(cmdstr);   
                
        outfile=sprintf('%s.predict_labels',predictdatafile);
        fid=fopen(outfile,'r');
        cycle=true;
        counter=0; 
        predictlabels=[];
        while cycle==true
            counter=counter+1;
            currentLine=fgets(fid);                
            if currentLine==-1%remove header of predict file
                break;
            end
            if counter>1
                t=str2num(currentLine);   
                if t(1)==1
                    predictlabels=[predictlabels;t(t(1)+1)];
                else
                    predictlabels=[predictlabels;-t(t(1)+1)];
               end
            end     
        end            
        fclose(fid);
                
        %output
        predict_values=[];
        predictinfofile=sprintf('%s%s.feature.vectors.txt.predict.info',params.relativePath,params.queryfile);   
        predictinfo=load(predictinfofile);
        index=0;
        for j=1:2:size(predictinfo,1)
            t=predictinfo(j,:);
            h=predictinfo(j+1,:);
            index=index+1;
            if t(2)==0 && h(2)==0
                predict_values=[predict_values;0];
            else                      
                if t(2)==1 && h(2)==1
                    if abs(predictlabels(j))>=abs(predictlabels(j+1))
                        predict_values=[predict_values;predictlabels(j)];
                    else
                        predict_values=[predict_values;predictlabels(j+1)];
                    end
                elseif t(2)==1
                        predict_values=[predict_values;predictlabels(j)];
                elseif h(2)==1
                        predict_values=[predict_values;predictlabels(j+1)];
                else
                        predict_values=[predict_values;0];
                end
            end
        end         
        finalpredictdatafile=sprintf('%s%s.predict.target.homolog.decvalues.txt',params.relativePath,params.queryfile); 
        dlmwrite(finalpredictdatafile,predictlabels,' ');
        finalpredictdatafile=sprintf('%s%s.predict.final.decvalues.txt',params.relativePath,params.queryfile); 
        dlmwrite(finalpredictdatafile,predict_values,' ');   
            
        
    
            
          
        
        
        
        
        