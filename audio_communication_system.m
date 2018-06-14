%%%%%%%%%%%%%%%%%%%
%%% Audio Input %%%
%%%%%%%%%%%%%%%%%%%

%wav = input('Plase input the filename of the audio: \n', 's');

% sample: datapoints of the audio
% fs: frequency of samplerate
% nbits: bits of sampling 
%[sample, fs, nbits] = wavread(wav);

% cnt_point: number of datapoints
% cnt_track: number of tunnels
% delta_t: interval of TD
% t: sequence of time for sampling
%[cnt_point, cnt_track] = size(sample);
%delta_t = 1 / fs;
%t = (0:1:cnt_point = 1) / fs;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate Random Sequence %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = 8;
K = log2(M);
sym_map=[1;(1+j)/sqrt(2);j;(-1+j)/sqrt(2);-1;(-1-j)/sqrt(2);-j;(1-j)/sqrt(2)]; %8PSK symbols

Ns = 10;              % Number of symbols
%bits = round(rand(K,Ns));           % KxNs matrix of random 0,1 bits
bits = [
     0     0     1     1     0     1     0     1     1     0;
     0     0     0     0     1     1     0     0     0     0;
     1     1     1     1     1     1     1     0     1     1];% For test

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate 8PSK Signal %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mpsk_sin, mpsk_cos] = mpsk(bits, Ns);
subplot(2,2,3)
plot(mpsk_sin);
subplot(2,2,4)
plot(mpsk_cos);
mpsk_signal = mpsk_sin + mpsk_cos;
figure(2)
plot(mpsk_signal);






Es = 10.^([[-7] [8:1:22]]/10); % Energy per symbol
Eb = Es/K;                % Energy per bit
na = length(Es);          % number of energy per symbol
No = 2;                   % noise unit variance (watt/Hz)
Es_No = Es/No;              % EsNo
Eb_No = Eb/No;              % EbNo

BER = zeros(1,na);
SER = zeros(1,na);
Pseint = zeros(1,na);

dphi = 0.01*pi/M;                   % interval of $\phi$
phi = [-pi/M+dphi/2:dphi:pi/M];     % $\phi$
nphi = length(phi);                 % number of $\phi$

