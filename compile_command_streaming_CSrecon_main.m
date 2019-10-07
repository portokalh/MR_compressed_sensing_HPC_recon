%compile me
function compile_command_streaming_CSrecon_main
%exec_env_var='';%optional shell env var to be cleared.
include_files = {'/cm/shared/workstation_code_dev/recon/CS_v2/gui_info_collect.m' ...    
    '/cm/shared/workstation_code_dev/recon/CS_v2/puller_glusterspaceCS_2.m' ...
    '/cm/shared/apps/MATLAB/R2015b/toolbox/signal/signal/hamming.m' ...
    '/cm/shared/apps/MATLAB/R2015b/toolbox/images/images/padarray.m' ...
    '/cm/shared/workstation_code_dev/recon/CS_v2/zpad.m'...
    };%optional, but required if using exec_env_var, can be empty.
compile_dir=compile_command__allpurpose('streaming_CS_recon_main_exec.m',include_files);%,exec_env_var);

code_dir=fileparts(mfilename('fullpath'));
original_builtin_script = 'run_streaming_CS_recon_main_exec_builtin_path.sh';
original_builtin_path=fullfile(code_dir,original_builtin_script);

% This has been modified to use a shell script which checks for CS_CODE_DEV
% env var, and runs the desired veresion. Defaults to stable if
% unspecified, when developing code, use latest.
% update_bin_cmd=sprintf('cp %s %s/;rm %s;ln -s %s/%s %s',original_builtin_path,compile_dir,bin_path,compile_dir,original_builtin_script,bin_path)
% system(update_bin_cmd);
update_bin_cmd=sprintf('cp %s %s/',original_builtin_path,compile_dir);
system(update_bin_cmd);
exit