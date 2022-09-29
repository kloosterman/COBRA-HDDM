%% quantile prob plots of observed vs simulated data

% simdata = readtable('/Users/kloosterman/Dropbox/PROJECTS/COBRA/hddm/123back_bias_novelvsfam/data/data_simulated.csv');
% obsdata = readtable('/Users/kloosterman/Dropbox/PROJECTS/COBRA/hddm/123back_bias_novelvsfam/data/data_observed.csv');

simdata = readtable('/Users/kloosterman/Dropbox/PROJECTS/COBRA/hddm/123back_bias_novelvsfam/data/data_simulated_basicmodel.csv');
obsdata = readtable('/Users/kloosterman/Dropbox/PROJECTS/COBRA/hddm/123back_bias_novelvsfam/data/data_observed_basicmodel.csv');

%% make quantiles observed data
% quantiles = [0 20 40 60 80 100];
quantiles = [10 20 40 60 80 90];
% quantiles = [10 30 50 70 90];
% quantiles = [0.5 10 30 50 70 90 99.5];
nbins = length(quantiles)-1;

SUBJ = unique(obsdata.subj_idx);

% for observed data
% includes stim non-target or target: want plots for stim present and
% absent
bindat = [];
for isub = 1:length(SUBJ)
  for istim = 0:1 % non-target or target
    for icond = 2:3 % only 2 and 3 back
      line_inds = obsdata.subj_idx == SUBJ(isub) & obsdata.stim == icond & obsdata.stimulus == istim;
      dat = [abs(obsdata.rt(line_inds)) obsdata.accuracy(line_inds) ]; % we want choice (0 is no), not accuracy
      prctile_edges = prctile(dat(:,1), quantiles);       % get RT quantiles
      [~,~,bin] = histcounts(dat(:,1), prctile_edges); % get bin that each trial belongs to
      dat = [dat bin]; % append      
      for ibin = 1:nbins % get proportion responses per quantile
        for icorr = 0:1 % incorrect and correct
          bindat(isub, istim+1, ibin, icond-1, icorr+1, 1 ) = mean(dat(dat(:,3)==ibin & dat(:,2) == icorr, 1)); %mean RT per bin and correct/incorrect
          bindat(isub, istim+1, ibin, icond-1, icorr+1, 2 ) = length(dat(dat(:,3)==ibin & dat(:,2) == icorr, 2)) ./ length(find(dat(:,3)==ibin)); %correctness
        end
      end
    end
  end
end

% for simulated data sets (50)
% includes stim non-target or target: want plots for stim present and
% absent
datasets_sim = unique(simdata.sample);

bindat_sim = [];
for idata = 1:10%length(datasets_sim)
  disp(idata)
  curdat = simdata(simdata.sample==datasets_sim(idata),:);
  for isub = 1:length(SUBJ)
    for istim = 0:1 % non-target or target
      for icond = 2:3 % only 2 and 3 back
        line_inds = strcmp( sprintf('wfpt(%d.%d).%d', icond, istim, SUBJ(isub)), simdata.node ) & simdata.sample == datasets_sim(idata); 
%         line_inds = simdata.subj_idx == SUBJ(isub) & simdata.stim == icond & simdata.stimulus == istim;
        dat = [abs(simdata.rt(line_inds )) simdata.response(line_inds) ]; 
        prctile_edges = prctile(dat(:,1), quantiles);       % get RT quantiles
        [~,~,bin] = histcounts(dat(:,1), prctile_edges); % get bin that each trial belongs to
        dat = [dat bin]; % append
        for ibin = 1:nbins % get proportion responses per quantile
          for icorr = 0:1 % incorrect and correct
            bindat_sim(isub, idata, istim+1, ibin, icond-1, icorr+1, 1 ) = median(dat(dat(:,3)==ibin & dat(:,2) == icorr, 1)); %mean RT per bin and correct/incorrect
            bindat_sim(isub, idata, istim+1, ibin, icond-1, icorr+1, 2 ) = length(dat(dat(:,3)==ibin & dat(:,2) == icorr, 2)) ./ length(find(dat(:,3)==ibin)); %correctness
          end
        end
      end
    end
  end
end

  
%% plot quantile prob plot
bindat_avg = squeeze(nanmean(bindat)); % stim bins cond incorrect/correct RT/accuracy
bindat_sim_avg = squeeze(nanmean(bindat_sim)); % stim bins cond incorrect/correct RT/accuracy
titnames = {'Target absent' 'Target present'};
SAV=1;
YLIM = [0.5 1.5];
f = figure; f.Position = [        1000        1115         625         150];
quantile_leg = {};
for istim = 1:2 % absent, present
  subplot(1,2,istim); hold on
  title(titnames{istim})
  xlabel('Proportion of responses');
  ylabel('Reaction time (s)')
  ylim(YLIM)
  % plot the observed data
  for ibin = 1:nbins
    
    dat_x = squeeze(bindat_avg(istim, ibin,:,:,2 )); % RT for both conditions, correct and incorrect
    dat_y = squeeze(bindat_avg(istim, ibin,:,:,1 )); % RT for both conditions, correct and incorrect
    dat_x(:,2) = flipud(dat_x(:,2)); % to get the line straight through the points
    dat_y(:,2) = flipud(dat_y(:,2));

    plot(dat_x(:), dat_y(:), '-x')
%     quantile_leg{ibin} = sprintf('q %1.2f', quantiles(ibin+1)/100-0.1);
    quantile_leg{ibin} = sprintf('q %d-%d', quantiles(ibin), quantiles(ibin+1) );
  end
  for idata = 1:size(bindat_sim_avg,1)
    % plot the sim data
    for ibin = 1:nbins
      
      dat_x = squeeze(bindat_sim_avg(idata,istim, ibin,:,:,2 )); % RT for both conditions, correct and incorrect
      dat_y = squeeze(bindat_sim_avg(idata,istim, ibin,:,:,1 )); % RT for both conditions, correct and incorrect
      dat_x(:,2) = flipud(dat_x(:,2)); % to get the line straight through the points
      dat_y(:,2) = flipud(dat_y(:,2));
      
      plot(dat_x(:), dat_y(:), '.k')
%       quantile_leg{ibin} = sprintf('q %1.2f', quantiles(ibin+1)/100-0.1);
    end
  end
  %   if istim==2
  legend(quantile_leg, 'Location', 'EastOutside'); legend boxon
  %   end
  ax=gca;
  ax.FontSize = 8;
%   line()
end

% legend off
if SAV
  PREOUT = '/Users/kloosterman/Dropbox/PROJECTS/COBRA/hddm/123back_bias_novelvsfam/data/plots_quantile_prob';
  saveas(f, fullfile(PREOUT, sprintf('Quantprob_%dbins.pdf', nbins )))
  saveas(f, fullfile(PREOUT, sprintf('Quantprob_%dbins.png', nbins )))
end

