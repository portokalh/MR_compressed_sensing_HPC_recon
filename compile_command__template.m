function compile_command__template
exec_env_var='';%optional shell env var to be cleared.
include_files = {};%optional, but required if using exec_env_var, can be empty.
compile_command__allpurpose('mfile',include_files,exec_env_var);
