function T = proc_timeresolvedDDMpars()

PREIN = '/Users/kloosterman/Dropbox/PROJECTS/COBRA/hddm/123back_bias_novelvsfam/data/';
cd(PREIN)

timenames = {'early' 'middle' 'late'};
parnames = {'a', 'v', 't', 'z', 'dc' };

pars=NaN(156,3,6,5); % dimord subj_time_nback_param
parcondnames = {};
for itime = 1:3
  csv_in = fullfile(PREIN, sprintf('params__biasmodel_%s', timenames{itime}));
  T = readtable(csv_in);
  for ipar = 1:5
    for iback=1:3
      pars(:,itime,iback,ipar) = T.mean(contains(T.Var1, sprintf('%s_subj(%d)', parnames{ipar}, iback))); % nsub
      parcondnames{itime,iback,ipar} = sprintf('%s_%s_%dback', parnames{ipar}, timenames{itime}, iback);
    end
    %3-2back
    pars(:,itime,4,ipar) = pars(:,itime,3,ipar) - pars(:,itime,2,ipar);
    parcondnames{itime,4,ipar} = [parcondnames{itime,3,ipar} '-2back'];
    %2-1back
    pars(:,itime,5,ipar) = pars(:,itime,2,ipar) - pars(:,itime,1,ipar);
    parcondnames{itime,5,ipar} = [parcondnames{itime,2,ipar} '-1back'];
    %3-1back
    pars(:,itime,6,ipar) = pars(:,itime,3,ipar) - pars(:,itime,1,ipar);
    parcondnames{itime,6,ipar} = [parcondnames{itime,3,ipar} '-1back'];
  end
  clear T
end
% TODO contrasts

% flatten and save to csv
pars = reshape(pars, 156, []);
parcondnames = reshape(parcondnames, 1, []);

T = array2table(pars, 'VariableNames', parcondnames);

outfile = 'COBRA_biasDDM_time_by_nback.csv';
writetable(T, fullfile(fileparts(PREIN), outfile))
fprintf('%s written to\n%s', outfile, fileparts(PREIN))


