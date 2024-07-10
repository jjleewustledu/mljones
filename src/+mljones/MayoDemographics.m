classdef MayoDemographics < handle
    %% line1
    %  line2
    %  
    %  Created 16-Apr-2024 23:42:38 by jjlee in repository /Users/jjlee/MATLAB-Drive/mljones/src/+mljones.
    %  Developed on Matlab 24.1.0.2537033 (R2024a) for MACA64.  Copyright 2024 John J. Lee.
    
    properties
        days_separation_tol
        workdir
    end

    methods
        function this = MayoDemographics(varargin)
            this.workdir = fullfile(getenv("MATLABDRIVE"), "mljones", "data");
            this.days_separation_tol = 365;
        end

        
        function T = age_info(this)
            T_ = this.table_mayo();
            sub = T_.PatientIdentifier_SequentialForMayoAndLONIPTIDForADNI_;
            age = T_.AgeAtScanTime_Years_;
            T = natsortrows(table(sub, age), [], 1);
        end

        function T = apoe4_info(this)
            T_ = this.table_mayo();
            sub = T_.PatientIdentifier_SequentialForMayoAndLONIPTIDForADNI_;
            allele1 = T_.APOEAllele1;
            allele2 = T_.APOEAllele2;
            apoe4 = (allele1 == 4) + (allele2 == 4);
            T = natsortrows(table(sub, apoe4));
        end

        function T = cdr_info(this)
            T_ = this.table_mayo();
            sub = T_.PatientIdentifier_SequentialForMayoAndLONIPTIDForADNI_;
            cdr = T_.ClinicalDementiaRatingScaleGlobalScore;
            T = natsortrows(table(sub, cdr), [], 1);
        end
        
        function T = globbed_info(this)
            ld = load(fullfile(this.workdir, "globbed423.mat"));
            globbed = ascol(ld.globbed423);
            sub = NaN(length(globbed), 1);
            for idx = 1:length(globbed)
                re = regexp(mybasename(globbed(idx)), "sub-(?<sub>\d{3})_\S+", "names");
                sub(idx) = str2double(re.sub);
            end
            Filelist = globbed;
            T = natsortrows(table(Filelist, globbed, sub), [], [2, 3]);
        end

        function T = lewy_info(this)
            T_ = this.table_mayo();
            sub = T_.PatientIdentifier_SequentialForMayoAndLONIPTIDForADNI_;
            lewy = T_.LewyPathology_0_None_1_Brainstem_2_Transitional_3_Diffuse_;
            T = natsortrows(table(sub, lewy), [], 1);
        end

        function T = metaroi_info(this)
            T_ = this.table_mayo();
            sub = T_.PatientIdentifier_SequentialForMayoAndLONIPTIDForADNI_;
            metaroi = T_.FDG_PETInJagust_ADNIROIs;
            T = natsortrows(table(sub, metaroi), [], 1);
        end

        function T = parkinson_info(this)
            T_ = this.table_mayo();
            sub = T_.PatientIdentifier_SequentialForMayoAndLONIPTIDForADNI_;
            parkinson = T_.UnifiedParkinsonDiseaseRatingScale;
            T = natsortrows(table(sub, parkinson), [], 1);
        end

        function T = remsbd_info(this)
            T_ = this.table_mayo();
            sub = T_.PatientIdentifier_SequentialForMayoAndLONIPTIDForADNI_;
            remsbd = T_.REMSleepBehaviorDisorder;
            T = natsortrows(table(sub, remsbd), [], 1);
        end
        
        function T = sex_info(this) 
            T_ = this.table_mayo();
            sub = T_.PatientIdentifier_SequentialForMayoAndLONIPTIDForADNI_;
            sex = cell(length(sub), 1);
            sex(logical(T_.Sex)) = {'M'};
            sex(~logical(T_.Sex)) = {'F'};
            % sex = ascol(categorical(sex));
            T = natsortrows(table(sub, sex), [], 1);
        end

        function T = table_esm(~)
            %% Cohort_Mayo_ADNI_ indicates:
            %  423 from Mayo, 410 from ADNI.

            T = readtable( ...
                fullfile(getenv("SINGULARITY_HOME"), "MAYO", "table_esm.csv"));
        end

        function T = table_fdg(this)
            %% globbed_info with matched cdr

            T_ = this.globbed_info;
            Filelist = T_.Filelist;
            globbed = T_.globbed;
            sub = T_.sub;
            cdr = this.cdr_info().cdr;
            age = this.age_info().age;
            sex = this.sex_info().sex;
            apoe4 = this.apoe4_info().apoe4;
            lewy = this.lewy_info().lewy;
            parkinson = this.parkinson_info().parkinson;
            remsbd = this.remsbd_info().remsbd;
            Metaroi = this.metaroi_info().metaroi;

            T = table(Filelist, globbed, sub, cdr, age, sex, apoe4, lewy, parkinson, remsbd, Metaroi); % age, sex, apoe4
        end

        function T = table_mayo(~)
            %% Cohort_Mayo_ADNI_ indicates:
            %  423 from Mayo, 410 from ADNI.

            T = readtable( ...
                fullfile(getenv("SINGULARITY_HOME"), "MAYO", "table_mayo.csv"));
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
