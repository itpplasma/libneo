function run_condor_neo2(kpath)
%##############################################################
%function run_condor_neo2(kpath)
%##############################################################
% description:
%--------------------------------------------------------------
% This scripts runs NEO-2 with condor.
%##############################################################
% input:
%--------------------------------------------------------------
% kpath      ... path of run
%##############################################################  

%author:   Philipp Ulbl
%created:  21.01.2020

        mpath = pwd();
        cd(kpath);
        
        start_time = datetime;
        disp(['Start of NEO-2 at ', datestr(start_time)]); %write start time
        [~, ~] = system('./run_condor.sh');

        %loop until finished
        while(true)
            pause(60); % check each minute

            %get status of condor
            [~, res] = system('condor_q');
            %extract number of jobs
            openjobs = extractBetween(res, 'Total for all users: ', ' jobs');
            runjobs = extractBetween(res, 'Total for all users: ', ' suspended');
            runjobs = extractBetween(runjobs, 'idle, ', ' running');
            totaljobs = sum(n_points);
            finishedjobs = totaljobs - str2double(openjobs{:});
            %disp status
            disp(['Time elapsed: ', char(datetime - start_time), ...
                  '; Jobs: ', num2str(finishedjobs), '/', num2str(totaljobs), ' finished; ', ...
                  runjobs, ' running.'])

            if(openjobs == 0)
                break;
            end
        end
        disp(['Finished NEO-2 at ', datestr(datetime)]); %write end time
        
        cd(mpath)
end