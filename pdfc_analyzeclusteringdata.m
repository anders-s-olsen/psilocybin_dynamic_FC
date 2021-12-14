function p_cell = pdfc_analyzeclusteringdata(resultscell,survtbl,options)
% Statistical models of clustering data including fractional occurrence and
% dwell time. This script contains call functions for R scripts for
% survival analysis and permutation testing in linear mixed-effects models.
%
% Input:
%%% resultscell - Cell array with information on each scan session and
% centroid, including fractional occurrences and covariate levels
%%% survtbl - table with information on all state occurrences including the
% number of active samples in row.
%%% options - struct with many fields as specified in pdfc_check_input.m
% Output:
%%% p_cell - cell array containing evaluated statistical models for every
%%% computed centroid

resultscell_hdr = resultscell(1,:);
covnames        = options.covnames;

% initialize p_cell
p_cell_hdr{1,1} = 'N_centroids';
p_cell_hdr{1,2} = 'centroid_cur';
for cov = 1:options.numcovs
    numcols = size(p_cell_hdr,2);
    p_cell_hdr{1,numcols+1}  = ['FO_pval_',covnames{cov}];
    p_cell_hdr{1,numcols+2}  = ['FO_pval_perm',covnames{cov}];
    p_cell_hdr{1,numcols+3}  = ['FO_estimate_',covnames{cov}];
    p_cell_hdr{1,numcols+4}  = ['FO_intercept_',covnames{cov}];
    p_cell_hdr{1,numcols+5}  = ['FO_CIlower_',covnames{cov}];
    p_cell_hdr{1,numcols+6}  = ['FO_CIupper_',covnames{cov}];
    p_cell_hdr{1,numcols+7}  = ['DT_pval_',covnames{cov}];
    p_cell_hdr{1,numcols+8}  = ['DT_estimate_',covnames{cov}];
    p_cell_hdr{1,numcols+9}  = ['DT_CIlower_',covnames{cov}];
    p_cell_hdr{1,numcols+10} = ['DT_CIupper_',covnames{cov}];
    p_cell_hdr{1,numcols+11} = ['DT_curves_',covnames{cov}];
end

p_cell  = p_cell_hdr;
counter = 2;
for cov = 1:options.numcovs
    
    % Run permutation testing in R
    if options.run_perm_FO(cov)
        
        clearvars tblperm
        tblperm = cell2table(resultscell(2:end,[...
            find(strcmp(resultscell_hdr,'Subject')),...
            find(strcmp(resultscell_hdr,'Session')),...
            find(strcmp(resultscell_hdr,'N_centroids')),...
            find(strcmp(resultscell_hdr,'Current_centroid')),...
            find(strcmp(resultscell_hdr,'Fractional_occupancy')),...
            find(strcmp(resultscell_hdr,covnames{cov}))]...
            ),'VariableNames',{'Subject','Session','N_centroids',...
            'Current_centroid','Frac_occ','cov1'});
        
        tblperm(isnan(tblperm.cov1),:) = [];
        writetable(tblperm,[options.functionpath,'perm_in.csv']);
        system(['Rscript ',options.functionpath,'pdfc_permmaxT.R']);
        perm_lme_occ = readtable([options.functionpath,'perm_out.csv']);
        
        delete([options.functionpath,'perm_in.csv'])
        delete([options.functionpath,'perm_out.csv'])
        
    end
    
    
    for k = options.min_k:options.max_k
        
        disp(['Running stats: k = ', num2str(k)])
        
        % find cell indices for current N_centroids
        i_N_centroids = strcmp(resultscell(:,strcmp(resultscell_hdr,'N_centroids')),num2str(k));
        
        if options.run_dwell_time
            surv_N_centroids = survtbl.N_centroids==k;
        end
        
        for centroid_cur = 1:k
            
            % find cell indices for current centroid
            i_centroid_cur = strcmp(resultscell(:,strcmp(resultscell_hdr,'Current_centroid')),num2str(centroid_cur));
            
            % find cell indices for current N_centroids and centroid
            i_samples_cur  = i_N_centroids&i_centroid_cur;
            
            if options.run_frac_occ(cov)
                subject        = resultscell(i_samples_cur,strcmp(resultscell_hdr,'Subject'));
                cov1           = [resultscell{i_samples_cur,strcmp(resultscell_hdr,covnames{cov})}]';
                frac_occupancy = [resultscell{i_samples_cur,strcmp(resultscell_hdr,'Fractional_occupancy')}]';
                
                tbl            = table(subject,cov1,frac_occupancy);
                tbl(isnan(tbl.cov1),:) = [];
                formula_occ    = 'frac_occupancy ~ cov1 + (1|subject)';
                
                % statistics frac_occ
                lme_occ        = fitlme(tbl,formula_occ);
            end
            if options.run_perm_FO(cov)
                
                kcur           = perm_lme_occ.k==k;
                statecur       = perm_lme_occ.state==centroid_cur;
                lme_occ_p_perm = perm_lme_occ.pval_perm(kcur&statecur);
                
            end
            
            
            % statistics dwell time
            if options.run_dwell_time
                surv_centroid_cur     = survtbl.start_state == centroid_cur;
                surv_tbl_cur          = survtbl(surv_N_centroids&surv_centroid_cur,:);
                
                surv_tbl_cur.cov1     = surv_tbl_cur.(covnames{cov});
                writetable(cell2table(num2cell(options.DTintervals(cov,:))'),[options.functionpath,'DTintervals.csv']);
                
                writetable(surv_tbl_cur,[options.functionpath,'dwell_time.csv']);
                system(['Rscript ',options.functionpath,'pdfc_cox_frailty.R']);
                dwell_time_stats      = readtable([options.functionpath,'dwell_time_Surv_stats.csv']);
                dwell_time_curves     = readtable([options.functionpath,'dwell_time_Surv_curves.csv']);
                
                % ensure that the wrong ones aren't read; delete files
                delete([options.functionpath,'dwell_time_Surv_stats.csv']);
                delete([options.functionpath,'dwell_time_Surv_curves.csv']);
                delete([options.functionpath,'DTintervals.csv']);
                delete([options.functionpath,'dwell_time.csv']);
                
                dt_curves             = [];
                for DTint             = 1:length(options.DTintervals(cov,:))
                    dt_curves(:,DTint)    = dwell_time_curves.(['Val_',num2str(DTint)]);
                end
                
                dwell_time_p          = dwell_time_stats.p;
            end
            
            
            p_cell{counter,strcmp(p_cell_hdr,'N_centroids')}                   = k; %N_centroids
            p_cell{counter,strcmp(p_cell_hdr,'centroid_cur')}                  = centroid_cur; % centroid_cur
            if options.run_frac_occ(cov)
                if options.run_perm_FO(cov)
                    p_cell{counter,strcmp(p_cell_hdr,['FO_pval_perm',covnames{cov}])}  = lme_occ_p_perm;
                end
                p_cell{counter,strcmp(p_cell_hdr,['FO_pval_',covnames{cov}])}      = lme_occ.Coefficients.pValue(2);
                
                p_cell{counter,strcmp(p_cell_hdr,['FO_estimate_',covnames{cov}])}  = lme_occ.Coefficients.Estimate(2);
                p_cell{counter,strcmp(p_cell_hdr,['FO_intercept_',covnames{cov}])} = lme_occ.Coefficients.Estimate(1);
                coefCI                                                             = lme_occ.coefCI;
                p_cell{counter,strcmp(p_cell_hdr,['FO_CIlower_',covnames{cov}])}   = coefCI(2,1);
                p_cell{counter,strcmp(p_cell_hdr,['FO_CIupper_',covnames{cov}])}   = coefCI(2,2);
            end
            if options.run_dwell_time
                p_cell{counter,strcmp(p_cell_hdr,['DT_pval_',covnames{cov}])}      = dwell_time_p;
                p_cell{counter,strcmp(p_cell_hdr,['DT_estimate_',covnames{cov}])}  = dwell_time_stats.beta;
                p_cell{counter,strcmp(p_cell_hdr,['DT_CIlower_',covnames{cov}])}   = dwell_time_stats.betaCI1;
                p_cell{counter,strcmp(p_cell_hdr,['DT_CIupper_',covnames{cov}])}   = dwell_time_stats.betaCI2;
                p_cell{counter,strcmp(p_cell_hdr,['DT_curves_',covnames{cov}])}    = dt_curves;
            end
            
            counter = counter + 1;
            
        end
        disp(['Done with stats for model k= ',num2str(k)])
    end
    
    counter = 2;
end
fprintf('Done with all stats!\n')