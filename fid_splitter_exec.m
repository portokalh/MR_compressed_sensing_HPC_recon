function [ jobids ] = fid_splitter_exec(local_fid )
% For GRE/mGRE CS scans, will carve up fid into one fid for each independent
% volume.
%   
% Written by BJ Anderson, CIVM
% 19 September 2017


local_fid = 
[ready,bhdr]=check_subvolume_ready_in_fid(input_fid,n_volumes,bbytes,scanner,user,options)




end

