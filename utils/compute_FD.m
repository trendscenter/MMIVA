function [FD_abs, FD_rms ] = compute_FD(mp,flag)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [FD_abs, FD_rms, FD_rmsvd ] = compute_FD(mp,afniflag)
%  mp: nT x 6 rigid body head movement estimates
% flag: 2 if motion parameters are estimated using FSL 
%       1  if motion parameters are estimated using AFNI  
%       0 otherwise (assumes estimated using SPM by default)     
%
% NOTE THAT THIS IS ONLY APPROPRIATE FOR MOTION ESTIMATES FROM SPM and AFNI
% convert rotations into translations
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    flag = 0;
end
rad = 50; % assume head radius of 50mm
if ischar(mp)
    mp = load(mp);
	mp = mp(:,1:6);
else
    mp = mp(:,1:6);
end


if flag == 1 % AFNI
    nT = size(mp,1);
    mp = mp(:,[5 6 4 2 3 1]).*[-ones(nT,1) -ones(nT,1) ones(nT,1) (pi/180)*ones(nT,1) -(pi/180)*ones(nT,1) -(pi/180)*ones(nT,1)];
elseif flag == 2 % FSL 
    mp = mp(:,[4:6 1:3]);
end

rot = mp(:,4:6);
rdis = rad*tan(rot);
mp(:,4:6) = rdis;

FD_abs = sum(abs(diff(mp)),2); % Power et al 2012
FD_rms = sqrt(sum((diff(mp)).^2,2)); % Eswar
% FD_rmsvd = sqrt(sum(mp(2:end,1:3).^2,2)) - sqrt(sum(mp(1:end-1,1:3).^2,2)); % Van Dijk et al 2012
