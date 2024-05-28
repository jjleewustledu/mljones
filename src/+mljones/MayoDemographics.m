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

        function T = table_esm(~)
            %% Cohort_Mayo_ADNI_ indicates:
            %  423 from Mayo, 410 from ADNI.

            T = readtable( ...
                fullfile(getenv("SINGULARITY_HOME"), "MAYO", "table_esm.csv"));
        end

        function T = table_fdg(this)
            %% globbed_info with matched cdr

            T_ = this.globbed_info;
            globbed = T_.globbed;
            Filelist = T_.Filelist;
            sub = T_.sub;
            ses = T_.ses;
            Tcdr = this.cdr_info();
            Tamy = this.amyloid_info();
            Tage = this.age_info();
            Tsex = this.sex_info();
            Tapoe4 = this.apoe4_info();

            warning("off", "MATLAB:badsubscript")
            cdr = NaN(size(sub));
            amyloidosis = NaN(size(sub));
            age = NaN(size(sub));
            sex = cell(size(sub));
            apoe4 = NaN(size(sub));
            for row = 1:size(T_,1)
                the_sub = sub(row);
                the_ses = ses(row);

                try % cdr
                    Ucdr = Tcdr(Tcdr.sub == the_sub, :);
                    dday = abs(Ucdr.ses - the_ses);
                    Ucdr = addvars(Ucdr, dday, NewVariableNames="dday");
                    Ucdr = sortrows(Ucdr, "dday");
                    if Ucdr.dday(1) < this.days_separation_tol
                        cdr(row) = Ucdr.cdr(1);
                    else
                        fprintf("sub %g, ses %g, cdr %g had delta days ~ %g\n", ...
                            the_sub, the_ses, Ucdr.cdr(1), Ucdr.dday(1))
                    end
                catch ME
                    handwarning(ME)
                end

                try % age
                    Uage = Tage(Tage.sub == the_sub, :);
                    dday = the_ses - Uage.ses;
                    Uage = addvars(Uage, dday, NewVariableNames="dday");
                    abs_dday = abs(Uage.ses - the_ses);
                    Uage = addvars(Uage, abs_dday, NewVariableNames="abs_dday");
                    Uage = sortrows(Uage, "abs_dday");
                    age(row) = Uage.age(1) + years(days(Uage.dday(1)));
                catch ME
                    handwarning(ME)
                end

                try % sex
                    Usex = Tsex(Tsex.unique_subs == the_sub, :);
                    sex{row} = Usex.sex(1);
                catch ME
                    handwarning(ME)
                end

                try % apoe4
                    Uapoe4 = Tapoe4(Tapoe4.unique_subs == the_sub, :);
                    apoe4(row) = Uapoe4.apoe4(1);
                catch ME
                    handwarning(ME)
                end

                try % amyloid, must be last task block in for-loop
                    Uamy = Tamy(Tamy.sub == the_sub, :);
                    if ~isempty(Uamy)
                        Uamy__ = Uamy; % chronological
                        dday = abs(Uamy.ses - the_ses);
                        Uamy = addvars(Uamy, dday, NewVariableNames="dday");
                        Uamy = sortrows(Uamy, "dday"); % ordered by abs delta day from the_ses
                        if Uamy.dday(1) < this.days_separation_tol % select any amloid_s1 within datetime_separation_tol of acqdate
                            amyloidosis(row) = Uamy.amyloidosis(1);
                            continue
                        end
                        if all(~Uamy__.amyloidosis) % no known amyloidosis through end of amy scanning
                            if the_ses <= max(Uamy__.ses) + this.days_separation_tol
                                amyloidosis(row) = false;
                                continue
                            else
                                continue
                            end
                        end
                        if Uamy__.amyloidosis(1) % known amyloidosis since start of amy scanning
                            if the_ses >= min(Uamy__.ses) - this.days_separation_tol
                                amyloidosis(row) = true;
                                continue
                            else
                                continue
                            end
                        end

                        % at least one amy+ scan after first amy scan
                        [~,idx] = max(Uamy__.amyloidosis);
                        amyloidosis(row) = the_ses >= Uamy__.ses(idx) - this.days_separation_tol;
                    end
                catch ME
                    handwarning(ME)
                end
            end
            warning("on", "MATLAB:badsubscript")

            %sex = categorical(sex);

            T = table(Filelist, globbed, sub, ses, cdr, amyloidosis, age, sex, apoe4); % age, sex, apoe4
            T = T(~isnan(cdr) & ~isnan(amyloidosis), :);
            T.amyloidosis = logical(T.amyloidosis);
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
