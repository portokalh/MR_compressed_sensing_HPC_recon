%'heike','N54818','N54817_01','ser10'

scanner='kamy';
user='omega';
study = 'S67650_04';
series = 'ser04';
runno = 'S67650_70';

workdir = ['/glusterspace/' runno '.work/'];

if ~exist(workdir, 'dir')
    system(['mkdir -m 777 ' workdir]);
end

datapath=['/home/mrraw/' study '/' series '.fid'];

current_fid_size=get_remote_filesize(datapath);

if (current_fid_size == 0)
    find_file_cmd=['ssh ' user '@' scanner ' "ls -tr /home/vnmr1/vnmrsys/exp*/acqfil/fid | tail -n1"'];
    [junk,result]=system(find_file_cmd);
    latest_fid=result(1:end-1);
    input_fid=latest_fid;
else
    input_fid=datapath;
end

local_hdr = [workdir runno '_hdr.fid'];

if ~exist(local_hdr,'file');
    get_hdr_from_fid(input_fid,local_hdr,scanner,user);
end

[npoints,nblocks,ntraces,bitdepth,bbytes,complete_file_size,~]
[~,nblocks,~,~,bbytes,final_file_size]=load_fid_hdrCS(local_hdr);

for vv = 1:nblocks
    required_file_size = 32+vv*bbytes
    current_fid_size=get_remote_filesize(input_fid,scanner)
    [ready,bhdr]=check_subvolume_ready_in_fid(input_fid,vv,bbytes,scanner,user)
    local_volume_fid = sprintf('%s%s%s%02i%s',workdir, runno, '_m', (vv-1), '.fid');
    if ~exist(local_volume_fid,'file')
        tt=0;
        while (~ready)
            [ready,bhdr]=check_subvolume_ready_in_fid(input_fid,vv,bbytes,scanner,user)
            current_cycle = tt
            pause(30)
            current_fid_size=get_remote_filesize(input_fid,scanner)
            tt=tt+1;
        end
        sprintf('%s%02i%s','fid now contains data for volume ', (vv-1), '; attempting to pull from scanner')
        get_subvolume_from_fid(input_fid,local_volume_fid,vv,bbytes,scanner,user);
    end
end