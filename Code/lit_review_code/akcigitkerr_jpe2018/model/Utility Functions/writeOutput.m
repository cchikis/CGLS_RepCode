function [] = writeOutput(param,paramNames,pAux,pANames,m,filename)

if isstruct(param)
   pS = cell2mat(struct2cell(param));    % this part used for estimation
else
   pS = param;
end

if isstruct(pAux)
   pA = cell2mat(struct2cell(pAux));
else
   pA = pAux;
end

% Printing the final results
% construct cell array
result      = cell(m.nMoment,6);
result(:,1) = num2cell(m.momentVec);
result(:,2) = num2cell(m.momentData);
result(:,3) = num2cell(1:m.nMoment);
result(:,4) = m.momentName;
result(:,5) = num2cell(m.momentWgt);
result(:,6) = num2cell(m.mErr);

 
fd = fopen(filename,'w+');
fprintf(fd,'%10s\n','ESTIMATION RESULTS');
fprintf(fd,'%10s %10s %5s %30s %10s %10s\n','model','data','#','description','weight','error');

for i=1:m.nMoment
   fprintf(fd,'%10.6f %10.6f %5i %30s %10.6f %10.6f\n',result{i,1},result{i,2},result{i,3},result{i,4},result{i,5},result{i,6});
end
fprintf(fd,'\n\n');
fprintf(fd,'score = %10.12f\n\n',m.score);
fprintf(fd,'%10s\n','PARAMETERS');
  
for i = 1: length(pS)
   fprintf(fd,'%8s : ',char(paramNames{i}));
   fprintf(fd,'%1.12f\n',pS(i));
end
fprintf(fd,'\n\n');
fprintf(fd,'%10s\n','Equilibrium Objects');
for i = 1: length(pA)
   fprintf(fd,'%8s : ',char(pANames{i}));
   fprintf(fd,'%1.12f\n',pA(i));
end
  

fprintf(fd,'\n\n');
fclose(fd);
  
end