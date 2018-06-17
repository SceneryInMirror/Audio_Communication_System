%% Audio Input

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

[sample, fs] = audioread('music.wav', [1 294828], 'native');
[cnt_point, cnt_track] = size(sample);
sample = sample(:, 1);
bits = reshape(sample', [3, cnt_point/3]);

M = 8;
K = log2(M);
sym_map=[1;(1+1i)/sqrt(2);1i;(-1+1i)/sqrt(2);-1;(-1-1i)/sqrt(2);-1i;(1-1i)/sqrt(2)]; %8PSK symbols

%% Generate Random Sequence

Ns = cnt_point/3;              % Number of symbols, strongly influence the CPU, max 10000
%bits = [
%     0     0     1     1     0     1     0     1     1     0;
%     0     0     0     0     1     1     0     0     0     0;
%     1     1     1     1     1     1     1     0     1     1];% For test
%bits = round(rand(K,Ns));           % KxNs matrix of random 0,1 bits
 
%% Generate 8PSK Signal

% WAY 1:
Nb = 10; % point number of carrier for one period, strongly influence the scatterplot
fc = 1; % frequency of carrier, also the frequency of symbol
t = 0:1/(Nb * fc):1/fc - 1/(Nb * fc); % time sequence for a period
carrier = exp(1i * 2 * pi * fc * t);
  
s_mpsk = [];
s_I = [];
s_Q = [];
test = [];

for n=1:Ns
    k = 4 * bits(1, n) + 2 * bits(2, n) + bits(3, n) + 1;
    test = [test k];
    s_mpsk = [s_mpsk real(sym_map(k) * carrier)]; % generate 8PSK signal
    s_I = [s_I real(sym_map(k)) * cos(2 * pi * fc * t)]; % the I branch of 8PSK
    s_Q = [s_Q imag(sym_map(k)) * sin(2 * pi * fc * t)]; % the Q branch of 8PSK
end

figure(3)
subplot(4,1,1)
plot(s_mpsk)
subplot(4,1,2)
plot(abs(fft(s_mpsk)))
subplot(4,1,3)
plot(s_I)
subplot(4,1,4)
plot(s_Q)

% WAY 2:
%[s_Q, s_I] = mpsk(bits, Ns);
%subplot(2,2,3)
%plot(s_Q);
%subplot(2,2,4)
%plot(s_I);
%mpsk_signal = s_I - s_Q;
%figure(2)
%subplot(2,1,1)
%plot(mpsk_signal);
%subplot(2,1,2)
%plot(abs(fft(mpsk_signal)));

%% Upsample

s_upsample = upsample(s_mpsk, 8); % 8 times upsample

figure(4)
subplot(2,1,1)
plot(s_upsample);
subplot(2,1,2)
plot(abs(fft(s_upsample)));

%% Lowpass Filter

%s_transmit = lowpass_transmit(s_upsample);
s_transmit = s_upsample; % for test

%% AWGN Channel

%Es = 10.^([[-7] [8:1:22]]/10); % Energy per symbol
Es = 10; %for test
Eb = Es/K;                % Energy per bit
N0 = 2;
SNR = 10 * log10((K * Eb/N0) / (8 * fc));
s_awgn = awgn(s_transmit, SNR, 'measured');

figure(5)
subplot(2,1,1)
plot(s_awgn);
subplot(2,1,2)
plot(abs(fft(s_awgn)));

%% Lowpass Filter

%s_receive = lowpass_transmit(s_awgn);
s_receive = s_awgn; %for test

%% Downsample

s_downsample = downsample(s_receive, 8); % 8 times downsample


%% 8PSK Judgement

s_demodulate_I = s_downsample .* cos(2 * pi * fc * repmat(t,[1, Ns]));
s_demodulate_Q = s_downsample .* sin(2 * pi * fc * repmat(t,[1, Ns]));

s_demodulate_I = reshape(s_demodulate_I, [Nb, Ns]);
s_demodulate_Q = reshape(s_demodulate_Q, [Nb, Ns]);

s_demodulate_I = 2.0 / Nb * sum(s_demodulate_I);
s_demodulate_Q = 2.0 / Nb * sum(s_demodulate_Q);

figure(6)
subplot(2,1,1)
plot(s_demodulate_I);
subplot(2,1,2)
plot(s_demodulate_Q);

figure(7)
plot(s_demodulate_I, s_demodulate_Q, 'b.')

s_result = s_demodulate_I - 1i * s_demodulate_Q;
distance = abs(repmat(s_result, [M, 1]) - repmat(sym_map, [1, Ns]));

[min_dis, min_pos] = min(uint32(distance .* 10000));

figure(8)
subplot(3,1,1)
plot(test)
subplot(3,1,2)
plot(min_pos)

subplot(3,1,3)
plot(test-min_pos)

min_pos = min_pos - 1;
bits_result = [];
bits_result = [bits_result sign(bitand(min_pos, 4))];
bits_result = [bits_result; sign(bitand(min_pos, 2))];
bits_result = [bits_result; mod(min_pos, 2)];

%% Calculate



%%
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


