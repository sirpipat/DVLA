function minfo = mooringinfo(name, metadatadir)
% minfo = MOORINGINFO(name, metadatadir)
%
% Get the info of any given CAATEX mooring
%
% INPUT:
% name          name of a mooring
% metadatadir   location where the metadata files are stored
%
% OUTPUT:
% minfo         info of that mooring
%
% EXAMPLE:
% info = mooringinfo('SIO3');
%
% SEE ALSO:
% LISTCAATEXMOORINGS
%
% Last modified by spipatprathanporn@ucsd.edu, 02/09/2026

% TODO: change to environment variables
% TODO: decide where to host data+metadata
defval('metadatadir', '/Volumes/aog/AOGquake/work/Ivy/CAATEX/metadata');

% determine which directories / metadata files to read
if strcmpi(name, 'SIO1')
    metadatafile = fullfile(metadatadir, 'minfo', 'minfo.SN101.mat');
    load(metadatafile, 'minfo');
    metadatafile = fullfile(metadatadir, 'minfo', 'minfo.SN102.mat');
    minfo2 = load(metadatafile, 'minfo');
elseif strcmpi(name, 'SIO2')
    metadatafile = fullfile(metadatadir, 'minfo', 'minfo.SN103.mat');
    load(metadatafile, 'minfo');
    metadatafile = fullfile(metadatadir, 'minfo', 'minfo.SN104.mat');
    load(metadatafile, 'minfo2');
elseif strcmpi(name, 'SIO3')
    metadatafile = fullfile(metadatadir, 'minfo', 'minfo.SN105.mat');
    load(metadatafile, 'minfo');
    metadatafile = fullfile(metadatadir, 'minfo', 'minfo.SN106.mat');
    load(metadatafile, 'minfo2');
elseif strcmpi(name, 'NERSC1')
    metadatafile = fullfile(metadatadir, 'minfo', 'minfo.SN110.mat');
    load(metadatafile, 'minfo');
elseif strcmpi(name, 'NERSC2')
    metadatafile = fullfile(metadatadir, 'minfo', 'minfo.SN131.mat');
    load(metadatafile, 'minfo');
elseif strcmpi(name, 'NERSC3')
    metadatafile = fullfile(metadatadir, 'minfo', 'minfo.SN133.mat');
    load(metadatafile, 'minfo');
else
    error('Invalid station name')
end

% merge info from 2 sets of hydrophones for SIO moorings
if contains(name, 'SIO')
    minfo.STARtag = {minfo.STARtag, minfo2.minfo.STARtag};
    minfo.dstardep = [minfo.dstardep minfo2.minfo.dstardep];
    minfo.xdcrdep = [minfo.xdcrdep minfo2.minfo.xdcrdep];
    minfo.snx = horzcat(minfo.snx, minfo2.minfo.snx);
    minfo.hmdep = horzcat(minfo.hmdep, minfo2.minfo.hmdep);
    minfo.indep = horzcat(minfo.indep, minfo2.minfo.indep);

% make STARtag a cell array for NERSC moorings to be consistent with SIO
else
    minfo.STARtag = {minfo.STARtag};
end
end