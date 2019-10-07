%compile me
function compile_command_setup
%exec_env_var='';%optional shell env var to be cleared.
include_files = {'/cm/shared/apps/MATLAB/R2015b/toolbox/stats/stats/quantile.m' ...
    '/cm/shared/workstation_code_dev/recon/CS_v2/CS_utilities/read_header_of_CStmp_file.m' ... 
    '/cm/shared/workstation_code_dev/recon/CS_v2/sparseMRI_v0.2/init.m'
    };%optional, but required if using exec_env_var, can be empty.
compile_command__allpurpose('setup_volume_work_for_CSrecon_exec.m',include_files);%,exec_env_var);
