function [ jobids,msg1,msg2 ] = dispatch_slurm_jobs( batch_file,slurm_options_string,optional_dependencies,optional_dependency_type )
% Custom handling of using the sbatch command
%
jobids=0;
dependencies='';
msg1='';
msg2='';

if exist(batch_file,'file') % batch_file could just be a naked command
    bf_array = strsplit(batch_file,'/');
    batch_file_name=bf_array(length(bf_array));
    batch_file_name=batch_file_name{1};
    bf_array(length(bf_array)) =[];
    default_dir = strjoin(bf_array,'/');
    if ~strcmp(default_dir(end),'/')
       default_dir = [default_dir '/']; 
    end
end


if exist('optional_dependencies','var') && ~isempty(optional_dependencies)
    
    od_array_1 = strsplit(strjoin(strsplit(optional_dependencies,','),':'),':');
    
    dep_test = str2double(od_array_1{1});
    
    if ~isnumeric(dep_test) || isnan(dep_test)
        optional_dependency_type = od_array_1{1};
        od_array_1(1)=[];
    end
    
    
    list='';
    if ~isempty(od_array_1)
         list =[':' strjoin(od_array_1,':')];
    end
    
    if ~exist('optional_dependency_type','var')
        optional_dependency_type = 'afterok';
    end
    optional_dependencies = list;

    dependencies = ['--dependency=' optional_dependency_type optional_dependencies];
elseif exist('optional_dependency_type','var' )
    dependencies = ['--dependency=' optional_dependency_type];
end

if exist('slurm_options_string','var')
    out_test_array = strsplit(slurm_options_string,'--out');
    out_test = 2 - length(out_test_array);
    
    if out_test % no --out specified, use directory of batch file
        if exist(batch_file,'file') % batch_file could just be a naked command
            %bf_array = strsplit(batch_file,'/');
            %bf_array(length(bf_array)) = [];
            %out_dir = strjoin(bf_array,'/');
            slurm_out_string = 'slurm-%j.out'; % trying to avoid premature interpretation of '%j'...
            out_string = sprintf(' --out=%s%s ',default_dir,slurm_out_string);
            slurm_options_string=[slurm_options_string out_string];
        end
    end
else
    slurm_options_string='';
    if exist(batch_file,'file') % batch_file could just be a naked command
        %bf_array = strsplit(batch_file,'/');
        %bf_array(length(bf_array)) =[];
        %out_dir = strjoin(bf_array,'/');
        slurm_out_string = 'slurm-%j.out'; % trying to avoid premature interpretation of '%j'...
        out_string = sprintf(' --out=%s%s ',default_dir,slurm_out_string);
        slurm_options_string=[slurm_options_string out_string];
    end
    % Add other defaults
end

sbatch_cmd = sprintf('sbatch %s %s %s',slurm_options_string,dependencies,batch_file);%['sbatch --requeue --mem=' mem ' -s -p ' queue ' ' slurm_options ' ' setup_dependency ' --job-name=' job_name ' --out=' batch_folder 'slurm-%j.out ' batch_file];
[~,msg]=system(sbatch_cmd);
msg_string = strsplit(msg,' ');
jobid = strtrim(msg_string{end});
msg1=msg;
%disp(msg)


if ~isempty(str2num(jobid)) % On successful schedule, schedule a backup job.
    jobids=jobid;
    
    %% Code for creating backup jobs in case originals fail. 25 May 2017, BJA
    dependencies = [' --dependency=afternotok:' jobid ];
    backup_sbatch_cmd =  sprintf('sbatch %s %s %s',slurm_options_string,dependencies,batch_file);
    [~,msg]=system(backup_sbatch_cmd);
    msg_string = strsplit(msg,' ');
    jobid_bu = strtrim(msg_string{end});
    msg2=msg;    

    if ~isempty(str2num(jobid_bu))
        jobids=[jobids ':' jobid_bu];
    end
    
    if exist(batch_file,'file')
        %% Rename batch files
        if  ~isempty(str2num(jobid))
            rename_sbatch_cmd = ['cp ' batch_file ' ' default_dir jobid '_' batch_file_name]; % changed 'mv' to 'cp'
            system(rename_sbatch_cmd);
            msg1='';
        end
        if   ~isempty(str2num(jobid_bu))
            rename_sbatch_cmd = ['mv ' batch_file ' ' default_dir jobid_bu '_backup_' batch_file_name]; % New code
            system(rename_sbatch_cmd);
            msg2='';
        end
    end
else
    if exist(batch_file,'file')
        if ~isempty(str2num(jobid))
            rename_sbatch_cmd = ['mv ' batch_file ' ' default_dir jobid '_' batch_file_name]; % changed 'mv' to 'cp'
            system(rename_sbatch_cmd);
            msg1='';
        end
    end
end

end

