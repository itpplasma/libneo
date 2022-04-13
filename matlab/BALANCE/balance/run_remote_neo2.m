function run_remote_neo2(configfile, kpath)
%##############################################################
%function run_remote_neo2(configfile, kpath)
%##############################################################
% description:
%--------------------------------------------------------------
% This scripts runs NEO-2 on computers specified in the config
% file.
%##############################################################
% input:
%--------------------------------------------------------------
% configfile ... path to config file
% kpath      ... path of run
%##############################################################

%author:   Philipp Ulbl
%created:  21.01.2020

    %minimum RAM per NEO2 run
    MINRAM = 30000;
    %minimum # of CPUs per NEO2 run
    MINCPU = 1;

    %fastest way: textscan
    
    fid = fopen(configfile);
    %scan string, int, int
    raw = textscan(fid, '%s %d %d');
    %get indices of comment lines: #
    ind = ~cellfun(@(c) strcmp(c(1), '#'), raw{1});
    %remove comment lines
    raw = cellfun(@(c) vertcat(c(ind)), raw, 'UniformOutput', false);
    %assign hostnames, cpus and ram
    host = raw{1};
    cpu = raw{2};
    mem = raw{3};
    
    %keep other solutions as comments because they are nice
    
%     %read_in much faster than importdata

%     raw = read_in(configfile);
%     raw = raw(~cellfun(@(c) strcmp(c(1), '#'), raw));
%     raw = cellfun(@(c) strsplit(c, ' '), raw, 'UniformOutput', false);
% 
%     host = cellfun(@(c) vertcat(c{1}), raw, 'UniformOutput', false);
%     cpu = cellfun(@(c) vertcat(str2double(c{2})), raw);
%     mem = cellfun(@(c) vertcat(str2double(c{3})), raw);
    
%     %import config file
%     conf = importdata(configfile);
%     %remove commented lines
%     rmind = cellfun(@(c) strcmp(c(1), '#'), conf.textdata);
%     conf.textdata(rmind) = [];
%     conf.data(rmind, :) = [];
%     %extract data
%     host = conf.textdata;
%     cpu = conf.data(:, 1);
%     mem = conf.data(:, 2);
    
    %calculate maximum jobs per host
    maxjob = min(floor(mem ./ MINRAM), floor(cpu ./ MINCPU));
    
    %construct queue of hosts to run
    queue = {};
    for k = 1:numel(maxjob)
        queue = [queue, repmat(host(k), 1, maxjob(k))];
    end
    
    mpath = pwd();
    cd(kpath);
    
    %get all jobs (sx.xxx folders in kpath)
    jobs = dir(kpath);
    jobs = jobs(3:end); %ignore . and .. in folder
    jobs = jobs([jobs.isdir]==true); %get directories
    jobs = jobs(cellfun(@(s) strcmp(s, 'TEMPLATE_DIR'), {jobs.name})==false); %ignore TEMPLATE_DIR
    jobs = {jobs.name}; %get only names
    totaljobs = numel(jobs); %count
    
    %initialization
    jobsrun = [];
    queuerun = [];
    jobsfinished = [];
    queuewait = queue;
    jobswait = jobs;
    
    start_time = datetime;
    disp(['Start of NEO-2 at ', datestr(start_time)]); %write start time
    
    %first check for already finished jobs
    finished = [];
    for k = 1:numel(jobswait)
        %if file done.out is created, job is done
        if(exist([jobswait{k}, '/done.out'], 'File') == 2)
            finished = [finished, k]; %add index to finished array
        end
    end
    jobsfinished = [jobsfinished, jobswait(finished)];
    jobswait(finished) = [];
    
    %loop until finished
    while(true)
        
        %if jobs run, check if finished
        if(numel(jobsrun) > 0)
            finished = [];
            %go through all jobs that run
            for k = 1:numel(jobsrun)
                %if file done.out is created, job is done
                if(exist([jobsrun{k}, '/done.out'], 'File') == 2)
                    finished = [finished, k]; %add index to finished array
                end
            end
            %if finished is empty, nothing happens
            %otherwise move hobs from run to finished and queue from run to
            %wait
            jobsfinished = [jobsfinished, jobsrun(finished)];
            queuewait = [queuewait, queuerun(finished)];
            jobsrun(finished) = [];
            queuerun(finished) = [];
        end

        %check if hosts are available and jobs wait to run
        if(numel(queuewait) > 0 && numel(jobswait) > 0)

            %start jobs
            start = [];
            for k = 1:min(numel(queuewait), numel(jobswait))
                %number of next job
                jobnum = numel(jobsfinished) + numel(jobsrun) + 1;
                disp(['Start job ', num2str(jobnum), '/', num2str(numel(jobs)), ': ', jobswait{k}, ' @', queuewait{k}])
                %start job over ssh on queue that is waiting
                %create log file with hostname and console output and
                %create done.out if job is finished
                system(['ssh ', queuewait{k}, ' ''cd ', kpath, jobswait{k}, ' ; hostname > log.txt ; ./neo_2.x >> log.txt 2>&1 ; touch done.out'' &']);
                %system(['ssh ', queuewait{k}, ' ''cd ', kpath, jobswait{k}, ' ; >done.out'' &']); %dummy for test
                %shift job to run and queue to run
                jobsrun = [jobsrun, jobswait(k)];
                queuerun = [queuerun, queuewait(k)];
                %save index of started job
                start = [start, k];
            end
            %remove job from wait and queue from wait
            queuewait(start) = [];
            jobswait(start) = [];

        end
            
        
        %disp status
        disp(['Time elapsed: ', char(datetime - start_time), ...
              '; Jobs: ', num2str(numel(jobsfinished)), '/', num2str(totaljobs), ' finished'])
        
        %if all jobs are done, break
        if(numel(jobsfinished) == totaljobs)
            break;
        end
        
        pause(30); % check each half minute
    end
    disp(['Finished NEO-2 at ', datestr(datetime)]); %write end time

    cd(mpath);
end