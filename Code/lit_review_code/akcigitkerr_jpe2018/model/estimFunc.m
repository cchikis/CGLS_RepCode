function [score,m,eq,qualityFirmA, nProdA, qual_StepA, stepA] = estimFunc(parin)
  
  global alg

  nan_score  = Inf;
  inf_score  = Inf;
  bnd_score  = Inf;
  fail_score = Inf;
 
  params = alg.pScale;
  params(alg.estParams) = parin.*alg.pScale(alg.estParams);
  alg.paramTried = params;
  
  alg.pub(1) = 1 - 1/params(2) + .0000001;   % sigma is between 0 and 1 - 1/psiTilde

  if (any(params < alg.plb) || any(params > alg.pub))
    fprintf(1,'BOUNDS ERROR\n');
    score = bnd_score;
  else
    p = putParamsStructure(params); 
    p = putParams(p);
    [~,~,eq,flag]= solveEq(p);

    if (flag == 0)
      fprintf(1,'EQ SOLVE FAILED\n');
      score = fail_score;
    else
      [score,m] = callMoment(p,eq);
      fprintf(1,'\nCurrent score =\t %f\n',score);
      
    end

    if (isnan(score))
      fprintf(1,'SCORE IS NAN\n');
      score = nan_score;
    end
    if (isinf(score))
      fprintf(1,'SCORE IS INF\n');
      score = inf_score;
    end

    if (score < alg.bestval)
      alg.lastparin = parin;
      alg.bestval   = score;
      fd = fopen(alg.summaryBestStep,'w+');
      fprintf(fd,'Best score =\t %f\n',alg.bestval);
      fprintf(fd,'%10s\n','Corresponding Parameters');
      fprintf(fd,'%16.12f\n,',params);
      fprintf(fd,'%10s\n',datestr(now));
      fclose(fd);
      save 'logs/lastStepMoment.mat' m;
    end
  end

  fprintf(1,'Best score =\t %f\n',alg.bestval);
end
