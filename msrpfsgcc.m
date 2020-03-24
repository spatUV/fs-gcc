function [FSGCCmat,lags,tpwin] = msrpfsgcc(xn,Nfft,B,M)
% MSRPFSGCC Computes the Frequency-Sliding Generalized-Cross Correlation
% matrix for the two-microphone signal XN
% Inputs:
%      XN:       Nx2 signal matrix, where N is the signal length.
%    NFFT:       FFT size.
%       B:       Spectral window size (in frequency bins).
%       M:       Spectral window hop (in frequency bins).
%
% Outputs:
%   FSGCCmat:    FS-GCC Matrix.
%       LAGS:    Delay corresponding to each row.
%      TPWIN:    Temporal equivalent of spectral window (used in WSVD-GCC)
%
% Copyright (C) 2020 Maximo Cobos


% Nyquist bin
Kmax = round(Nfft/2) + 1;

% Number of frequency bands
L = max(1,fix((Kmax-B+M)/M)+1);

% Create spectral window
wind = hann(B,'periodic');
spwin = zeros(Nfft,1);
spwin(1:B/2) = wind(B/2+1:end);
spwin(end-B/2+1:end) = wind(1:B/2);
% Inverse spectral window
tpwin = real(ifft(spwin));

% Compute cross-power spectrum
tf  = fft(xn,Nfft);
CPS = tf(:,1).*conj(tf(:,2));
PHAT = exp(1i*angle(CPS));

% Get sub-band GCCs
FSGCCmat = zeros(Nfft,L);
for k = 0:L-1
    PHATd = circshift(PHAT, [-k*M, 0]);
    PHATd = PHATd.*spwin;
    FSGCCmat(:,k+1) = fftshift(ifft(PHATd)); 
end

lags = -(Nfft/2):Nfft/2-1;