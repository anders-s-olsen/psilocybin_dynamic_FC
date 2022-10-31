function ptbl = pdfc_analyzeclusteringdata(resultstbl,survtbl,options)
% ptbl = pdfc_analyzeclusteringdata(resultstbl,survtbl,options)
% Statistical models of clustering data including fractional occurrence and
% % dwell time. This script contains call functions for R scripts for
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
%
% Anders S Olsen April - November 2021, October 2022
% Neurobiology Research Unit, Copenhagen University Hospital Rigshospitalet


ptbl_variablenames = {'N_centroids','centroid_cur'};
for cov = 1:options.numcovs
    ptbl_variablenames  = [ptbl_variablenames,['FO_pval_',       options.covnames{cov}]];
    ptbl_variablenames  = [ptbl_variablenames,['FO_pvalcor',   options.covnames{cov}]];
    ptbl_variablenames  = [ptbl_variablenames,['FO_estimate_',   options.covnames{cov}]];
    ptbl_variablenames  = [ptbl_variablenames,['FO_intercept_',  options.covnames{cov}]];
    ptbl_variablenames  = [ptbl_variablenames,['FO_CIlower_',    options.covnames{cov}]];
    ptbl_variablenames  = [ptbl_variablenames,['FO_CIupper_',    options.covnames{cov}]];
    ptbl_variablenames  = [ptbl_variablenames,['DT_pval_',       options.covnames{cov}]];
    ptbl_variablenames  = [ptbl_variablenames,['DT_estimate_',   options.covnames{cov}]];
    ptbl_variablenames  = [ptbl_variablenames,['DT_CIlower_',    options.covnames{cov}]];
    ptbl_variablenames  = [ptbl_variablenames,['DT_CIupper_',    options.covnames{cov}]];
end
for cov = 1:options.numcovs
    ptbl_variablenames  = [ptbl_variablenames,['DT_curves_',     options.covnames{cov}]];
end


ptbl_variableclass = [{'int16','int16'},repelem({'double'},numel(ptbl_variablenames)-2-options.numcovs),repelem({'cell'},options.numcovs)];

ptbl  = table('Size',[0,numel(ptbl_variablenames)],...
    'VariableNames', ptbl_variablenames,...
    'VariableTypes',ptbl_variableclass);

counter = 1;
for cov = 1:options.numcovs
    
    % Run permutation testing in R
    if options.run_perm_FO(cov)
        
        clearvars tblperm
        tblperm = resultstbl;
        tblperm.cov1 = resultstbl.(options.covnames{cov});
        
        tblperm(isnan(tblperm.cov1),:) = [];
        writetable(tblperm,[options.functionpath,'perm_in.csv']);
        system(['cd ',options.functionpath,'; Rscript ',options.functionpath,'pdfc_permmaxT.R']);
        perm_lme_occ = readtable([options.functionpath,'perm_out.csv']);
        
        delete([options.functionpath,'perm_in.csv'])
        delete([options.functionpath,'perm_out.csv'])
        
    end
    
    
    for k = options.min_k:options.max_k
        
        disp(['Running stats for ',options.covnames{cov},', k = ', num2str(k)])
        
        % find cell indices for current N_centroids
        i_N_centroids = resultstbl.N_centroids==k;
        
        if options.run_dwell_time
            surv_N_centroids = survtbl.N_centroids==k;
        end
        
        for centroid_cur = 1:k
            
            % find tbl indices for current centroid
            i_centroid_cur = resultstbl.Current_centroid==centroid_cur;
            
            % find tbl indices for current N_centroids and centroid
            i_samples_cur  = i_N_centroids&i_centroid_cur;
            
            if options.run_frac_occ(cov)
                subject        = resultstbl.Subject(i_samples_cur);
                cov1           = resultstbl.(options.covnames{cov})(i_samples_cur);
                frac_occupancy = resultstbl.Fractional_occupancy(i_samples_cur);
                
                stattbl            = table(subject,cov1,frac_occupancy);
                stattbl(isnan(stattbl.cov1),:) = [];
                formula_occ    = 'frac_occupancy ~ cov1 + (1|subject)';
                
                % statistics frac_occ
                lme_occ        = fitlme(stattbl,formula_occ);
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
                surv_tbl_cur.cov1     = surv_tbl_cur.(options.covnames{cov});
                
                writetable(cell2table(num2cell(options.DTintervals(cov,:))'),[options.functionpath,'DTintervals.csv']);
                writetable(surv_tbl_cur,[options.functionpath,'dwell_time.csv']);
                
                system(['cd ',options.functionpath,'; Rscript ',options.functionpath,'pdfc_cox_frailty.R']);
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
                
            end
            
            ptbl.N_centroids(counter)                   = k; %N_centroids
            ptbl.centroid_cur(counter)                  = centroid_cur; % centroid_cur
            if options.run_frac_occ(cov)
                if options.run_perm_FO(cov)
                    ptbl.(['FO_pvalcor',options.covnames{cov}])(counter)  = lme_occ_p_perm;
                end
                ptbl.(['FO_pval_',options.covnames{cov}])(counter)      = lme_occ.Coefficients.pValue(2);
                
                ptbl.(['FO_estimate_',options.covnames{cov}])(counter)  = lme_occ.Coefficients.Estimate(2);
                ptbl.(['FO_intercept_',options.covnames{cov}])(counter) = lme_occ.Coefficients.Estimate(1);
                coefCI = lme_occ.coefCI;
                ptbl.(['FO_CIlower_',options.covnames{cov}])(counter)   = coefCI(2,1);
                ptbl.(['FO_CIupper_',options.covnames{cov}])(counter)   = coefCI(2,2);
            end
            if options.run_dwell_time
                ptbl.(['DT_pval_',options.covnames{cov}])(counter)      = dwell_time_stats.p;
                ptbl.(['DT_estimate_',options.covnames{cov}])(counter)  = dwell_time_stats.beta;
                ptbl.(['DT_CIlower_',options.covnames{cov}])(counter)   = dwell_time_stats.betaCI1;
                ptbl.(['DT_CIupper_',options.covnames{cov}])(counter)   = dwell_time_stats.betaCI2;
                ptbl.(['DT_curves_',options.covnames{cov}]){counter}    = dt_curves;
            end
            
            counter = counter + 1;
            
        end
%         disp(['Done with stats for model k= ',num2str(k)])
    end
    
    counter = 1;
end
fprintf('Done with all stats!\n')