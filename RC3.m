function Hd = RC3
%RC3 Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.4 and DSP System Toolbox 9.6.
% Generated on: 19-Jun-2018 22:59:29

% FIR Window Raised-cosine filter designed using the FIRRCOS function.

% All frequency values are in Hz.
Fs = 682666.66667;  % Sampling Frequency

N    = 100;           % Order
Fc   = 42666.666667;  % Cutoff Frequency
TM   = 'Rolloff';     % Transition Mode
R    = 0.4;           % Rolloff
DT   = 'Normal';      % Design Type
Beta = 0.5;           % Window Parameter

% Create the window vector for the design algorithm.
win = kaiser(N+1, Beta);

% Calculate the coefficients using the FIR1 function.
b  = firrcos(N, Fc/(Fs/2), R, 2, TM, DT, [], win);
Hd = dfilt.dffir(b);

% [EOF]
