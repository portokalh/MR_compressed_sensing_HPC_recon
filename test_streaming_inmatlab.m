% list of tests run by james to check on option handling and things, 
% This file is disposeable, but facilitated testing the many things.
% a test which should fail becuase itnlim is done twice, once by position,
% then again by name.
%streaming_CS_recon_main_exec kamy S00001t S67930_03 ser12 CS256_8x_pa18_pb54 1 Itnlim=2

% just setting intlim by position
%streaming_CS_recon_main_exec kamy S00001t S67930_03 ser12  1
% just setting cs_table by position
%streaming_CS_recon_main_exec kamy S00001t S67930_03 ser12 CS256_8x_pa18_pb54 
% setting intlim and cs_table by position
% streaming_CS_recon_main_exec kamy S00001t S67930_03 ser12 CS256_8x_pa18_pb54 1
% setting intlim and cs_table by postion, but in wrong order. This works,
% and i've decided that is okay.
%streaming_CS_recon_main_exec kamy S00001t S67930_03 ser12 1 CS256_8x_pa18_pb54
% streaming_CS_recon_main_exec kamy S00001t S67930_03 ser12 CS_table=CS256_8x_pa18_pb54 Itnlim=1

% streaming_CS_recon_main_exec kamy S00003t S67930_03 ser12 CS_table=CS256_8x_pa18_pb54 Itnlim=1 fid_archive

% streaming_CS_recon_main_exec kamy S00003t S67930_03 ser12 CS_table=CS256_8x_pa18_pb54 Itnlim=1 first_volume=2 last_volume=2

%streaming_CS_recon_main_exec a b c d help
% streaming_CS_recon_main_exec kamy S00001t S67930_03 ser12 CS_table=CS256_8x_pa18_pb54 Itnlim=1  

%streaming_CS_recon_main_exec kamy S00003t S67930_03 ser12 CS_table=CS256_8x_pa18_pb54 Itnlim=1 fid_archive

% streaming_CS_recon_main_exec kamy S00003t S67930_03 ser12 CS_table=CS256_8x_pa18_pb54 Itnlim=1 fid_archive overwrite
