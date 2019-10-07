 function [] = CS_recon_cluster_bj_multithread_v2(reconfile,volume_number,outpath,scale_file,target_machine,fermi_filter,chunk_size,CS_recon_params,sampling_fraction)

% Fermi_filter functionality beyond binary off/on with defaults is
% currently not supported.  The code below is in development, and requires
% that the cleanup exec function be successfully recompiled with variable
% fermi filter parameter functionality. BJA, 8 Nov 2016

if ~exist('fermi_filter','var')
    fermi_filter=0;
else
    if ischar(fermi_filter)
        fermi_filter_opts = strsplit(fermi_filter,'_');
        num_opts = length(fermi_filter_opts);
        fermi_filter = str2double(fermi_filter_opts{1}); % str2num if str2double fails
        if num_opts > 1
            w1 = str2double(fermi_filter_opts{2});
            if num_opts > 2
                w2 = str2double(fermi_filter_opts{3});
            end
        end
    else % If fermi_filter is just a scalar value...
        if fermi_filter
            w1 = 0.15;
            w2 = 0.75;
        end
    end
end


if (exist('CS_recon_params','var') && ~ischar(CS_recon_params))
    fermi_filter_opts = strsplit(fermi_filter,'_');
end


if ~exist('matlab_path','var')
    matlab_path = '/cm/shared/apps/MATLAB/R2015b/'; 
end

if ~exist('queue','var')
    queue='matlab';
end

% Default values for CS L1 Recon
TVWeight = 0.0012; 	% Weight for TV penalty - only TV is on, but I encourage you to try wavelets as well.
xfmWeight = 0.006;	% Weight for Transform L1 penalty
Itnlim = 48;		% Number of iterations

if (exist('CS_recon_params','var') && ischar(CS_recon_params))
    CS_recon_params_2 = strjoin(strsplit(CS_recon_params,' '),'');
    CS_recon_opts = strsplit(CS_recon_params_2,',');
    for ii=1:length(CS_recon_opts)
        var_and_val = strsplit(CS_recon_opts{ii},'=');
        if length(var_and_val) == 2
            variable_i = var_and_val{1};
            value_i = var_and_val{2};
            value_ii = str2double(value_i);
            switch variable_i
                case 'TVWeight'
                    if ((value_ii > 0) && (value_ii <= 1))
                        TVWeight = value_ii;
                    end
                case 'xfmWeight'
                    if ((value_ii > 0) && (value_ii <= 1))
                        xfmWeight = value_ii;
                    end
                case 'Itnlim'
                    if (value_ii) % Zero iterations not allowed!
                        Itnlim = ceil(abs(value_ii));
                    end
                otherwise
            end
        end
    end
    msg = ['Custom CS recon parameter(s) specified: TVWeight=' num2str(TVWeight) ...
        ', xfmWeight=' num2str(xfmWeight) ', and Itnlim=' num2str(Itnlim) '.'];
    
    disp(msg);
    disp('Please check that these values are correct.');
    disp(['(Your input was: ''' CS_recon_params '''.']);
end

wavelet_dims=[12 12]; % Can vary if need be.

%% Load data

load(reconfile);

%% Immediately check to see if volume has already been reconstructed
if  ~exist( 'do_russway','var')
    myvolstr =sprintf(['%0' num2str(numel(num2str(nvols-1))) 'i' ],volume_number-1);
else
    myvolstr = num2str(volume_number-1);
    if volume_number-1 < 10
        myvolstr = ['00' myvolstr];
    elseif volume_number-1 < 100
        myvolstr = ['0' myvolstr];
    end
end

dir1 = [runno '_m' myvolstr];
voldir = fullfile(outpath,dir1,[dir1 'images']);
if ~exist(voldir,'dir')
    mkdir(voldir);
end


if exist('procpar_path','var')
    [~,pp_name,pp_ext] = fileparts(procpar_path);
    procpar_for_archive = [voldir '/' pp_name pp_ext];

    if ~exist(procpar_for_archive,'file')
        cp_cmd = ['cp ' procpar_path ' ' voldir '/'];
        system(cp_cmd);
     end
end

volume_message = ['Current volume = ' myvolstr '.'];
%disp(
volume_message%)

finished_slices = dir( [voldir '/*.raw' ]);
finished_slices_count = length(finished_slices(not([finished_slices.isdir])));

should_I_do_work=0;
if (finished_slices_count ~= voldims(3))
    should_I_do_work=1;
else
    msg = ['All slices previously reconstructed.  Skipping CS recon for ' voldir '.'];
    msg
    %disp(msg)
end

if should_I_do_work
    
    batch_folder = [outpath '/' dir1 '/sbatch/'];
    if ~exist(batch_folder,'dir')
        dir_cmd = ['mkdir ' batch_folder];
        system(dir_cmd);
    
        chmod_cmd = ['chmod 777 ' batch_folder];
        system(chmod_cmd)
    end
    
    %% Create temporary volume for intermediate work
    % First [dims(1)] bytes form a header indicating which slices have already
    % been reconstructed. Also, we assume that dims1(1) = dims(1) always.
    
    header_size = voldims(1);
    work_folder = [outpath '/' dir1 '/work/']; % Will need to delete after work is completed
    
    if ~exist(work_folder,'dir')
        dir_cmd2 = ['mkdir ' work_folder '; chmod 777 ' work_folder];
        system(dir_cmd2);
    else
        chmod_cmd = ['chmod 777 ' work_folder];
        system(chmod_cmd);
    end
    
    temp_file = [work_folder '/' dir1 '.tmp'];
    
    work_done=zeros([header_size 1]);
    
    if exist(temp_file,'file')
        fid=fopen(temp_file,'r');
        work_done=fread(fid,voldims(1),'*uint8');
        fclose(fid);
    end
    
    slices_to_process = find(~work_done)';
    
    num_s2p = length(slices_to_process);
    
    
    %msg =
    ['Reconstructing ' num2str(num_s2p) ' slices.']%;
    
    %disp(msg)
    
    she_bang = '#!/bin/bash';
    slurm_options = '';%' --propagate=STACK ';
    
    reservation = getenv('CS_reservation');
    if (reservation)
        slurm_options = [slurm_options ' --reservation=' reservation ' '];
    end
    
    voldir = fullfile(outpath,dir1,[dir1 'images']);
    headfile = [voldir '/' dir1 '.headfile'];
    
    %% Determine if we need to construct our work.mat file
    %(or reconstruct in the oft-encountered case of a corrupt work.mat file)
    volume_variable_file = [work_folder dir1 '_workspace.mat'];
    run_setup = 0;
    try
        dummy = load(volume_variable_file,'aux_param.maskSize'); % Need to try to load an arbitrary variable from the work file
    catch
        run_setup = 1;
        if exist(volume_variable_file,'file');
            cmd=['rm ' volume_variable_file]
            system(cmd);
        end
    end
    
    
    if run_setup
        %% Make variable file
        variables_file = [work_folder dir1 '_setup_variables.mat'];
        
        if ~exist(variables_file,'file')
            variables.reconfile = reconfile;
            variables.procparpath = procpar_path;
            variables.outpath = outpath;
            variables.scale_file = scale_file;
            variables.target_machine = target_machine;
            variables.wavelet_dims = wavelet_dims;
            %variables.wavelet_type = wavelet_type;
            variables.TVWeight = TVWeight;
            variables.xfmWeight=xfmWeight;
            variables.Itnlim = Itnlim;
            if exist('sampling_fraction','var')
                variables.sampling_fraction = sampling_fraction;
            end
            save(variables_file,'variables');
        end
        
        local_fidpath=[outpath '/' runno '_m' myvolstr '.fid']; %%NEW
        
        if 0 % ~exist(local_fidpath,'file')
           %find out how many bytes per point
           if strcmp(bitdepth,'int16');
               bytes_per_point=2;
           else
               bytes_per_point=4;
           end
           %{
           tmp_fid_hd=[work_folder '/' runno '_fid_header.fid'];
           %fseek to the right place, skip 60 byte header and 28 byte block header, then data
            hd_dd_cmd=[ 'dd if=' fidpath ' of=' local_fidpath ' skip=0 count=1  bs=60'];
            display(hd_dd_cmd);
            system(hd_dd_cmd);
           
            byteskip=60+max_blocks*npoints*ntraces*bytes_per_point*(volume_number-1)+28*(volume_number-1)*(max_blocks);
            byteread=npoints*ntraces*bytes_per_point;
            %}

            %byteskip=32+max_blocks*npoints*ntraces*bytes_per_point*(volume_number-1)+28*(volume_number-1)*(max_blocks);
            byteread=28+npoints*ntraces*bytes_per_point;
            
            
            %{
            div=gcd(byteskip,byteread); 
            dd_cmd= [ 'dd if=' fidpath ' of=' local_fidpath ' skip=' num2str(byteskip/div) ' count=' num2str(byteread/div) ' bs=' num2str(div) ' conv=notrunc oflag=append'];
            %}
            
            dd_cmd = ['( dd bs=32 count=1 of=' local_fidpath ' && dd bs=' num2str((volume_number-1)*byteread) ' count=1 of=dev_null && dd bs=' num2str(byteread)  ' count=1 of=' local_fidpath ' conv=notrunc oflag=append ) < ' fidpath ];
            
            display(dd_cmd);
            system(dd_cmd);
  
        end
        %% Schedule setup via slurm and record jobid for dependency scheduling.
        
        if ~exist('setup_mem','var')
            setup_mem='50000';%'8196'; % Hopefully 8 Gb of memory is good enough for all setup jobs...will need to verify.
        end
        
        if ~exist('setup_executable_path','var')
            %setup_executable_path = '/glusterspace/BJ/BK/run_CS_recon_cluster_setup_work_exec.sh'; % Last stable: AX
            setup_executable_path = '/cm/shared/workstation_code_dev/recon/CS_v2/testing_and_prototyping/CS_recon_setup_executable/20170908_0903/run_CS_recon_cluster_setup_work_scratch_exec.sh'; % Last stable: BK % trying to add MGRE support DP
        end
        
        
        batch_file_name = [dir1 '_setup_work_for_CS_recon.bash'];
        batch_file = [batch_folder batch_file_name];
        
        job_name=[dir1 '_CS_recon_setup'];
        
        fid=fopen(batch_file,'w');
        fprintf(fid,'%s\n',she_bang);
        
        main_cmd = [setup_executable_path ' ' matlab_path ' '  variables_file  ' ' num2str(volume_number) ];
        
        fprintf(fid,'%s;\n',main_cmd);
        fclose(fid);
        setup_queue='high_priority';
        sbatch_cmd = ['sbatch --mem=' setup_mem ' -s -p ' setup_queue slurm_options ' --job-name=' job_name ' --out=' batch_folder 'slurm-%j.out ' batch_file];
        [~,msg]=system(sbatch_cmd);
        msg_string = strsplit(msg,' ');
        setup_jobid = strtrim(msg_string{end});
        
        rename_sbatch_cmd = ['mv ' batch_file ' ' batch_folder setup_jobid '_' batch_file_name]; % changed 'mv' to 'cp'
        system(rename_sbatch_cmd);

    end
    %%
    %     mask = skipint2skiptable(procpar_path); %sampling mask
    %     % % padarray to 512*512 for better recon?
    %     dyadic_idx = 2.^(1:12); %dyadic_idx = 2.^[1:12];
    %     pidx = find(max(voldims(2:3))<=dyadic_idx,1);
    %     p = 2^pidx;
    %     mask = padarray(mask,[p-voldims(2) p-voldims(3)]/2,0,'both');
    %     dims1 = [voldims(1) size(mask)];
    %
    %
    %     % Define recon parameters
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     % L1 Recon Parameters
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %     N = dims1(2:3); 		% image Size
    %     DN = N;		 	% data Size
    %     TVWeight = 0.0012; 	% Weight for TV penalty - only TV is on, but I encourage you to try wavelets as well.
    %     xfmWeight = 0.006;	% Weight for Transform L1 penalty
    %     Itnlim = 48;	%10	% Number of iterations
    %     OuterIt = length(TVWeight);
    %
    %
    if ~exist('mem','var')
        mem='5900'; % Want to allow 32-40 jobs per node, but use --ntasks-per-core=1 to make sure that every core has exactly one job on them.
    end
    
    if ~exist('cleanup_executable_path','var')
        %cleanup_executable_path = '/glusterspace/BJ/run_CS_recon_cluster_bj_cleanup_exec.sh';
        cleanup_executable_path = '/glusterspace/BJ/BD/run_CS_recon_cluster_bj_cleanup_and_fermi_filter_exec.sh';
        %cleanup_executable_path = '/home/rmd22/Documents/MATLAB/MATLAB_scripts_rmd/CS/run_CS_recon_cluster_bj_cleanup_and_fermi_filter_exec.sh';
    end
    
    if ~exist('executable_path','var')
        %executable_path = '/glusterspace/BJ/run_CS_recon_cluster_bj_slice_exec_V5.sh'; %V5 %V3
        executable_path = '/glusterspace/BJ/DZ/run_CS_recon_cluster_bj_slice_exec_V5.sh';
    end
    
    
    zero_width = ceil(log10((voldims(1)+1)));
    
    setup_dependency ='';
    if (run_setup && (str2double(setup_jobid)>1))
        setup_dependency = [' --dependency=afterok:' setup_jobid ];
    end
    
    
    
    num_chunks = ceil(length(slices_to_process)/chunk_size);
    %disp(
    ['Number of chunks per volume = ' num2str(num_chunks) '.']%)
    
    new_size = num_chunks*chunk_size;
    temp_size=length(slices_to_process);
    while new_size > temp_size
        slices_to_process = [slices_to_process NaN];
        temp_size = size(slices_to_process);
    end
    
    slices_to_process = reshape(slices_to_process,[chunk_size num_chunks]);
    
    for slice = slices_to_process
        slice_string = sprintf(['' '%0' num2str(zero_width) '.' num2str(zero_width) 's'] ,num2str(slice(1)));
        slice(isnan(slice))=[];
        if length(slice)>3
            no_con_test = sum(diff(diff(slice)));
        else
            no_con_test = 1;
        end
        
        for ss = 2:length(slice)
            if (no_con_test)
                slice_string = [slice_string '_' sprintf(['' '%0' num2str(zero_width) '.' num2str(zero_width) 's'] ,num2str(slice(ss)))];
            elseif (ss==length(slice))
                slice_string = [slice_string '_to_' sprintf(['' '%0' num2str(zero_width) '.' num2str(zero_width) 's'] ,num2str(slice(ss)))];
            end
        end
        
        batch_file_name = [dir1 '_slice' slice_string '_CS_recon.bash'];
        batch_file = [batch_folder batch_file_name];
        
        job_name=[dir1 '_CS_recon'];
        
        fid=fopen(batch_file,'w');
        fprintf(fid,'%s\n',she_bang);
        
        main_cmd = [executable_path ' ' matlab_path ' '  volume_variable_file  ' ' slice_string ];
        
        fprintf(fid,'%s;\n',main_cmd);
        fclose(fid);
        
        if chunk_size > 1
            plural = 's';
        else
            plural = '';
        end
        
        job_name=[job_name '_' num2str(chunk_size) '_slice' plural '_per_job'];
        
        %sbatch_cmd = ['sbatch --mem=' mem ' --ntasks-per-core=1 -s -p ' queue ' ' slurm_options ' ' setup_dependency ' --job-name=' job_name ' --out=' batch_folder 'slurm-%j.out ' batch_file];
        sbatch_cmd = ['sbatch --mem=' mem ' -c 2 --hint=nomultithread -p ' queue ' ' slurm_options ' ' setup_dependency ' --job-name=' job_name ' --out=' batch_folder 'slurm-%j.out ' batch_file];
        %sbatch_cmd = ['sbatch --mem=' mem ' -s -p ' queue ' ' slurm_options ' ' setup_dependency ' --job-name=' job_name ' --out=' batch_folder 'slurm-%j.out ' batch_file];
      
        %disp(sbatch_cmd)
        %         [~,msg]=system(sbatch_cmd);
        %         jobid = msg((end-6):(end-1));
        %
        [~,msg]=system(sbatch_cmd);
        msg_string = strsplit(msg,' ');
        jobid = strtrim(msg_string{end});
        
        %% Code for creating backup jobs in case originals fail. 25 May 2017, BJA
        backup_dependency = [' --dependency=afternotok:' jobid ];
        sbatch_cmd = ['sbatch --requeue --mem=' mem ' -s -p ' queue ' ' slurm_options ' ' backup_dependency ' --job-name=' job_name ' --out=' batch_folder 'slurm-%j.out ' batch_file];
        [~,msg]=system(sbatch_cmd);
        msg_string = strsplit(msg,' ');
        jobid_bu = strtrim(msg_string{end});
        
        %% End new code
        
        rename_sbatch_cmd = ['cp ' batch_file ' ' batch_folder jobid '_' batch_file_name]; % changed 'mv' to 'cp'
        system(rename_sbatch_cmd);
        
        rename_sbatch_cmd = ['mv ' batch_file ' ' batch_folder jobid_bu '_backup_' batch_file_name]; % New code
        system(rename_sbatch_cmd);
    end
    
    time_to_dispatch_slice_jobs=toc
    %disp(['Time to dispatch slice jobs = ' time_to_dispatch_slice_jobs '.'])
    
    %aux_param2.scaling=scaling;[slice_string sprintf(['' '%0' num2str(zero_width) '.' num2str(zero_width) 's'] ,num2str(slice(:)));
    aux_param2.dims=voldims;
    aux_param2.voldir=voldir;
    aux_param2.tempFile=temp_file;
    aux_param2.outpath=outpath;
    aux_param2.dir1=dir1;
    aux_param2.targetMachine = target_machine;
    aux_param2.fermi_filter=fermi_filter;
    aux_param2.headfile=headfile;
    
    
    if exist('w1','var')
        aux_param2.fermi_filter_w1=w1;
    end
    if exist('w2','var')
        aux_param2.fermi_filter_w2=w2;
    end
    
    cleanup_variable_file = [work_folder dir1 '_cleanup_variable.mat'];
    if ~exist(cleanup_variable_file,'file')
        %save(cleanup_variable_file,'struct1');
        save(cleanup_variable_file,'aux_param2');%,'-append');
    end
    batch_file_name = [dir1 '_CS_recon_cleanup.bash'];
    batch_file = [batch_folder batch_file_name];
    if ~exist('job_name','var')
        job_name=[dir1 '_CS_recon_cleanup'];
    end
    
    fid=fopen(batch_file,'w');
    fprintf(fid,'%s\n',she_bang);
    
    main_cmd = [cleanup_executable_path ' ' matlab_path ' '  cleanup_variable_file ];
    
    fprintf(fid,'%s;\n',main_cmd);
    fclose(fid);
    
    dependency ='';
    slices_to_process(isnan(slices_to_process))=[];
    if slices_to_process
        %dependency = ' --dependency=singleton --qos=90000';
        dependency = ' --dependency=singleton ';
    end
    
    queue= ' high_priority ';
    mem2 = '66000'%'32768'; % 32 GB
    sbatch_cmd = ['sbatch --mem=' mem2 ' -s -p ' queue ' ' slurm_options ' ' dependency ' --job-name=' job_name ' --out=' batch_folder 'slurm-%j.out ' batch_file];
    
    % Please note that " -t 15 " option which sets a 15 minute time limit
    % on the cleanup job.  It is possible that we will reach array sizes in
    % which this (somewhat arbitrarily imposed) constaint might need to be
    % either lifted or modified. BJA 8 Nov 2016
    % Upped to 30 for 512x512x1024 arrays
    
    %     [~,msg]=system(sbatch_cmd);
    %     jobid = msg((end-6):(end-1));
    
    [~,msg]=system(sbatch_cmd);
    msg_string = strsplit(msg,' ');
    jobid = strtrim(msg_string{end});
    
    %rename_sbatch_cmd = ['mv ' batch_file ' ' batch_folder jobid '_' batch_file_name];
    %system(rename_sbatch_cmd);
    
    %% Code for creating backup cleanup jobs in case originals fail; also serves as jobs to watch for the singleton runno_cleanup job. 13 July 2017, BJA
    runno_cleanup_job_name = ['runno_cleanup_for_' runno ];
    cleanup_backup_dependency = [' --dependency=afternotok:' jobid ];
    sbatch_cmd = ['sbatch --requeue --mem=' mem '  -s -p ' queue ' ' slurm_options ' ' cleanup_backup_dependency ' --job-name=' runno_cleanup_job_name ' --out=' batch_folder 'slurm-%j.out ' batch_file];
    [~,msg]=system(sbatch_cmd);
    msg_string = strsplit(msg,' ');
    jobid_bu = strtrim(msg_string{end});
    
    %% End new code
    
    rename_sbatch_cmd = ['cp ' batch_file ' ' batch_folder jobid '_' batch_file_name]; % changed 'mv' to 'cp'
    system(rename_sbatch_cmd);
    
    rename_sbatch_cmd = ['mv ' batch_file ' ' batch_folder jobid_bu '_backup_' batch_file_name]; % New code
    system(rename_sbatch_cmd);
    
    
    disp('Done sending CS recon work to the cluster for the current volume.')
end
end
