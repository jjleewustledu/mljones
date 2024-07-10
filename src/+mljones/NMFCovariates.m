classdef NMFCovariates < mladni.NMFCovariates & handle
    %% line1
    %  line2
    %  
    %  Created 16-Apr-2024 22:51:06 by jjlee in repository /Users/jjlee/MATLAB-Drive/mljones/src/+mljones.
    %  Developed on Matlab 24.1.0.2537033 (R2024a) for MACA64.  Copyright 2024 John J. Lee.
    
    methods
        function this = NMFCovariates(varargin)
            this = this@mladni.NMFCovariates(varargin{:});

            this.data_home_ = fullfile(getenv("SINGULARITY_HOME"), "MAYO");
            this.demogr_ = mljones.MayoDemographics();
            this.pet_on_T1w_suffix_ = 'orient-rpi.nii.gz';
        end

        function t = apply_table_qc(~, t)
        end

        function t = table_covariates(this)
            %% generates table in this.covariates_file() ~ 
            %  "Singularity/MAYO/NMF_FDG/baseline_cn/NumBases24/components/NMFCovariates_table_covariates_longitudinal.mat",
            %  for use by patterns_of_neurodegeneration_20240618.Rmd
            %
            %  sum(tc.Cohort == "AD") ~ 231
            %  sum(tc.Cohort == "PD") ~ 117
            %  sum(tc.Cohort == "REMSBD") ~ 43
            %  sum(tc.Cohort == "DLB") ~ 32
            %
            %  sum(tc.parkinson > 0 ) ~ 165
            %  sum(tc.remsbd > 0 ) ~ 59
            %  sum(tc.lewy > 0 ) ~ 32
            %
            %  sum(tc.remsbd > 0 & tc.lewy > 0) ~ 16
            %  sum(tc.remsbd > 0 & tc.parkinson > 0) ~ 42
            %
            %  sum(tc.remsbd > 0 & tc.parkinson > 0 & tc.lewy > 0) ~ 16

            if ~isempty(this.table_covariates_cache_)
                fprintf("%s: using cached in memory\n", stackstr())
                t = this.table_covariates_cache_;
                return
            end
            cache_file = this.covariates_file();
            if isfile(cache_file)
                fprintf("%s: using cached from filesystem\n", stackstr())
                ld = load(cache_file);
                this.table_covariates_cache_ = ld.t;
                t = this.table_covariates_cache_;
                return
            end

            % from MayoDemographics ~ 423 rows
            t = table_fdg(this); 

            % Cohort ~ categorical
            Cohort = repmat({'unknown'}, size(t,1), 1);
            Cohort(t.cdr == 0.5 & isnan(t.parkinson)) = {'CDR=0.5'};
            Cohort(t.cdr == 0.5 & t.parkinson == 0) = {'CDR=0.5'};
            Cohort(t.cdr == 0.5 & t.parkinson > 0) = {'CDR=0.5,UPDRS>0'};
            Cohort(t.cdr > 0.5 & isnan(t.parkinson)) = {'CDR>0.5'};
            Cohort(t.cdr > 0.5 & t.parkinson == 0) = {'CDR>0.5'};
            Cohort(t.cdr > 0.5 & t.parkinson > 0) = {'CDR>0.5,UPDRS>0'};
            assert(~any(contains(string(Cohort), 'unknown')))
            Cohort = categorical(Cohort);
            
            t = addvars(t, Cohort, NewVariableNames={'Cohort'});

            % Components ~ 423 rows
            t_ = this.table_selectedComponentWeightedAverageNIFTI;
            t = this.addvars_by_filelist(t, t_, t_.Components, NewVariableNames={'Components'});
            t = splitvars(t, 'Components');

            % apply table qc
            t = this.apply_table_qc(t);

            % store cache, & save/write table
            this.table_covariates_cache_ = t;
            save(cache_file, 't');
            writetable(t, strrep(cache_file, ".mat", ".csv"));
        end

        function t = table_fdg(this)
            if ~isempty(this.table_fdg_)
                t = this.table_fdg_;
                return
            end
            this.table_fdg_ = this.demogr_.table_fdg();
            t = this.table_fdg_;
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
