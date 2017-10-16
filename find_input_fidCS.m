function [ input_fid, local_or_streaming_or_static ] = find_input_fidCS( scanner,runno,study,series,user )
%
% local_or_streaming_or_static: 1 => local fid found, 2 -> using fid of in-progress scan, 3 -> using fid in its remote static location;

if ~exist('user','var')
    user='omega';
end

local_or_streaming_or_static = 3;

scratch_drive = getenv('BIGGUS_DISKUS'); % Currently, this is /glusterspace.
workdir = fullfile(scratch_drive,[runno '.work/']);
local_fid = [workdir runno '.fid'];

if exist(local_fid,'file')
    input_fid = local_fid;
    local_or_streaming_or_static = 1;
else
    datapath=['/home/mrraw/' study '/' series '.fid/fid'];
    current_fid_size=get_remote_filesize(datapath,scanner);
    
    if (current_fid_size == 0)
        find_file_cmd=['ssh ' user '@' scanner ' "ls -tr /home/vnmr1/vnmrsys/exp*/acqfil/fid | tail -n1"'];
        [~,result]=system(find_file_cmd);
        latest_fid=result(1:end-1);
        input_fid=latest_fid;
        local_or_streaming_or_static = 2;
    else
        input_fid=datapath;
    end
    
end
    
end

