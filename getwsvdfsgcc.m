function [GCCwsvd,lags] = getwsvdfsgcc(FSGCCmat,tpwin,lagmax)
% GETWSVDFSGCC retrieves a denoised GCC based on (Weighted) SVD, WSVD
% from a FS-GCC matrix 
% Inputs:
%   FSGCCMAT:          FS-GCC Matrix as computed by MSRPFSGCC.
%   LAGMAX (Optional): Restrict FS-GCC matrix to maximum lag.
%    
% Outputs:
%   GCCWSVD:      Recovered WSVD FS-GCC.
%
% Copyright (C) 2020 Maximo Cobos

N = size(FSGCCmat,1);
tpwin = fftshift(tpwin);

if ~exist('lagmax','var')
    lags = -N/2:N/2-1;
else
    FSGCCmat = FSGCCmat(N/2-lagmax:N/2+lagmax,:);
    tpwin = tpwin(N/2-lagmax:N/2+lagmax);
    lags = -lagmax:lagmax;  
    N = length(lags);
end

% SVD Low-rank approximation

ra1 = mean(abs(tpwin));
ra0 = sqrt(pi/2)*sqrt((tpwin'*tpwin)/N);
w = (ra0 - mean(abs(FSGCCmat)))./(ra0-ra1);
w(w<0)=0;
W = diag(w);

wsvdinput = FSGCCmat*W;
[U,~,~] = svd(wsvdinput,'econ');
GCCwsvd = real(U(:,1));
[~,maxabs] = max(abs(GCCwsvd));
GCCwsvd = GCCwsvd*sign(GCCwsvd(maxabs));
%GCCwsvd = GCCwsvd./norm(GCCwsvd);

