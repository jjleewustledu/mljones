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

        function t = apply_table_qc(this, t)
        end

        function t = table_covariates(this, opts)
        end

        function t = table_covariates_1stscan(this)
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
