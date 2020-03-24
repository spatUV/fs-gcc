function [GCCsvd,lags] = getsvdfsgcc(FSGCCmat,lagmax)
% GETSVDFSGCC retrieves a denoised GCC based on SVD from a FS-GCC matrix
% Inputs:
%   FSGCCMAT:          FS-GCC Matrix as computed by MSRPFSGCC.
%   LAGMAX (Optional): Restrict FS-GCC matrix to maximum lag.
%    
% Outputs:
%   GCCSVD:      Recovered SVD FS-GCC.
%
% Copyright (C) 2020 Maximo Cobos

N = size(FSGCCmat,1);

if ~exist('lagmax','var')
    lags = -N/2:N/2-1;
else
    FSGCCmat = FSGCCmat(N/2-lagmax:N/2+lagmax,:);
    lags = -lagmax:lagmax;  
end

% SVD Low-rank approximation
[U,~,~] = svd(FSGCCmat(:,:,1),'econ');
GCCsvd = real(U(:,1));
[~,maxabs] = max(abs(GCCsvd));
GCCsvd = GCCsvd*sign(GCCsvd(maxabs));
%GCCwsvd = GCCwsvd./norm(GCCwsvd);
