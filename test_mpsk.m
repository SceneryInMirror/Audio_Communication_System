%% Audio Input

%clear
M = 8;
K = log2(M);
sym_map=[1;(1+1i)/sqrt(2);1i;(-1+1i)/sqrt(2);-1;(-1-1i)/sqrt(2);-1i;(1-1i)/sqrt(2)]; %8PSK symbols

% sample: datapoints of the audio
% fs: frequency of samplerate
% cnt_point: number of datapoints
% cnt_track: number of tunnels
basepath = '/home/ykx/Audio_Communication_System/';
wav = [basepath 'music.wav'];
[sample, fs] = audioread(wav, 'native');
info = audioinfo(wav);
[cnt_point, cnt_track] = size(sample);

%sample = sample(:,1); % only use mono-channel statics
sample = reshape(sample, [], 1);
%sample = repmat(sample, [2,1]);
bits = reshape(double(dectobin(sample, 16)) - 48, 3, []);
Ns = length(bits);
BER = [];   
theory = [];

% display & save
figure(1)
subplot(3,1,1)
plot(sample)
title('Audio\_TD');
subplot(3,1,2)
plot(abs(fft(sample)))
title('Audio\_FD');
subplot(3,1,3)
plot(reshape(bits, 1, []))
title('Generated Binary Sequence');
saveas(1, [basepath 'figure/audio.png']);



%% Generate Random Sequence

%Ns = 100;              % Number of symbols, strongly influence the CPU, max 10000
%bits = [
%     0     0     1     1     0     1     0     1     1     0;
%     0     0     0     0     1     1     0     0     0     0;
%     1     1     1     1     1     1     1     0     1     1];% For test
%bits = round(rand(K,Ns));           % KxNs matrix of random 0,1 bits


%% Generate 8PSK Signal

fc = 1; % frequency of carrier, also the frequency of symbol

s_mpsk = [];
s_I = [];
s_Q = [];

k = 4 * bits(1, :) + 2 * bits(2, :) + bits(3, :) + 1;
s_mpsk = [s_mpsk sym_map(k)];
s_I = [s_I real(sym_map(k))]; % the I branch of 8PSK
s_Q = [s_Q imag(sym_map(k))]; % the Q branch of 8PSK
s = [s_I s_Q];

% display and save
figure(2)
subplot(2,1,1)
stem(s_I, '.')
title('Signal\_I')
subplot(2,1,2)
stem(s_Q, '.')
title('Signal\_Q')
saveas(2, [basepath 'figure/8psk.png']);


%% Upsample

s_upsample = upsample(s, 8);

% display and save
figure(3)
subplot(4,1,1)
stem(s_upsample(:,1), '.');
title('Upsample\_I\_TD')
subplot(4,1,2)
stem(abs(fft(s_upsample(:,1))), '.');
title('Upsample\_I\_FD')
subplot(4,1,3)
stem(s_upsample(:,2), '.');
title('Upsample\_Q\_TD')
subplot(4,1,4)
stem(abs(fft(s_upsample(:,2))), '.');
title('Upsample\_Q\_FD')
saveas(3, [basepath 'figure/upsample.png']);


%% Lowpass Filter

hd = HD;
myfilter = hd.Numerator;
s_transmit = filter(myfilter, 1, [s_upsample; zeros((length(myfilter) - 1) / 2, 2)]);
s_transmit = s_transmit(((length(myfilter)-1)/2+1): end, 1:2);

% display and save
figure(4)
subplot(2,1,1)
plot(s_transmit);
legend('signal\_I', 'signal\_Q');
title('Shaping\_TD')
subplot(2,1,2)
%stem(abs(fft(s_transmit)), '.');
plot(abs(fft(s_transmit)));
legend('signal\_I', 'signal\_Q');
title('Shaping\_FD')
saveas(4, [basepath 'figure/shaping.png']);

% combine two branches
s_transmit = s_transmit(:,1) + 1i * s_transmit(:,2);


%% AWGN Channel

Es_range = 10.^([[-7] [8:1:22]]/10); % Energy per symbol
%Es = 10; %for test
for item = 1:16
    Es = Es_range(item);
    Eb = Es/K;                % Energy per bit
    N0 = 2;
    SNR = 10 * log10((K * Eb/N0) / 8);
    s_awgn = awgn(s_transmit, SNR, 'measured');
    
    % display and save
    figure(5)
    subplot(2,1,1)
    plot(abs(s_awgn), 'b');
    hold on;
    plot(abs(s_transmit), 'r');
    legend('after\_awgn', 'before\_awgn');
    title('AWGN\_Amplitude');
    hold off;
    subplot(2,1,2)
    plot(angle(s_awgn), 'b');
    hold on;
    plot(angle(s_transmit), 'r');
    legend('after\_awgn', 'before\_awgn');
    title('AWGN\_Angle');
    hold off;
    saveas(5, fullfile(basepath, 'figure/', string(item), '/awgn.png'));
    
    %% Lowpass Filter

    %s_receive = real(s_awgn);
    %s_receive = [s_receive imag(s_awgn)];
    %s_receive = s_awgn;
    %s_receive = filter(myfilter, 1, [s_receive; zeros((length(myfilter) - 1) / 2, 2)]);
    %s_receive = s_receive(((length(myfilter)-1)/2+1): end, 1:2);
    s_receive = filter(myfilter, 1, [s_awgn; zeros((length(myfilter) - 1) / 2, 1)]);
    s_receive = s_receive(((length(myfilter)-1)/2+1): end, 1);
    %s_receive = real(s_awgn);
    %s_receive = [s_receive imag(s_awgn)];
    
    figure(6)
    subplot(3,1,1)
    plot(abs(s_receive), 'b');
    title('after filter');
    subplot(3,1,2)
    plot(abs(s_awgn), 'r');
    title('after awgn');
    subplot(3,1,3)
    plot(abs(s_transmit), 'y');
    title('before awgn');
    saveas(6, fullfile(basepath, 'figure/', string(item), '/matchedfilter.png'));


    %% Downsample

    s_downsample = downsample(s_receive, 8);
    s_downsample = [real(s_downsample) imag(s_downsample)];
    
    figure(7)
    %subplot(2,1,1)
    periodogram(s_downsample);
    %stem(s_downsample, '.');
    title('Downsample\_pdf');
    %subplot(2,1,2)
    %stem(abs(fft(s_downsample)), '.');
    %title('Downsample\_FD');
    saveas(7, fullfile(basepath, 'figure/', string(item), '/downsample.png'));


    %% 8PSK Judgement

    %s_demodulate_I = s_downsample_I;
    %s_demodulate_Q = s_downsample_Q;
    s_demodulate = s_downsample;

    %figure(8)
    %subplot(2,1,1)
    %plot(s_demodulate);
    %subplot(2,1,2)
    %plot(abs(fft(s_demodulate)));

    figure(8)
    plot(s_demodulate(:,1), s_demodulate(:,2), 'b.')
    title('Constellation Diagram');
    saveas(8, fullfile(basepath, 'figure/', string(item), '/constellation.png'));

    s_result = s_demodulate(:,1) - 1i * s_demodulate(:,2);
    s_result = s_result';
    distance = abs(repmat(s_result, [M, 1]) - repmat(sym_map, [1, Ns]));

    [min_dis, min_pos] = min(uint32(distance .* 10000));

    figure(9)
    subplot(3,1,1)
    plot(k)
    subplot(3,1,2)
    plot(min_pos)
    subplot(3,1,3)
    plot(k-min_pos)
    saveas(9, fullfile(basepath, 'figure/', string(item), '/error.png'));

    min_pos = min_pos - 1;
    bits_result = [];
    bits_result = [bits_result sign(bitand(min_pos, 4))]    ;
    bits_result = [bits_result; sign(bitand(min_pos, 2))];
    bits_result = [bits_result; mod(min_pos, 2)];

    figure(11)
    subplot(2,1,1)
    plot(reshape(bits, 1, []), 'b');
    title('transmitted\_bits');
    subplot(2,1,2)
    plot(reshape(bits_result, 1, []), 'r');
    title('received\_bits');
    saveas(11, fullfile(basepath, 'figure/', string(item), '/bits.png'));
    
    
    %% Calculate
    bits_error = abs(bits_result-bits);
    %BER = [BER sum(sum(bits_error(:,:) == 1))];
    BER = [BER biterr(bits, bits_result)];
    %theory = [theory 1/2 * erfc(1/2*10.^(SNR/10))];
    
end

theory = erfc(sqrt(Es_range/N0 * (sin(pi/M))^2));
figure(10)
hold off;
BE = BER;
BER = BER / (3 * length(bits));
plot(10 * log10(Es_range/(K * N0)) , 10 * log10(BER), 'b');
hold on;
%plot(10 * log10(Es_range/(K * N0)) , 10 * log10(theory), 'r');
plot(ebno0, 10 * log10(ber0), 'r');
title('BER(dB)-EbN0(dB)')
legend('simulate result', 'theoretical calculation');
saveas(10, [basepath 'figure/BER-EBN0.png']);

%%
%na = length(Es);          % number of energy per symbol
%No = 2;                   % noise unit variance (watt/Hz)
%Es_No = Es/No;              % EsNo
%Eb_No = Eb/No;              % EbNo

%BER = zeros(1,na);
%SER = zeros(1,na);
%Pseint = zeros(1,na);

%dphi = 0.01*pi/M;                   % interval of $\phi$
%phi = [-pi/M+dphi/2:dphi:pi/M];     % $\phi$
%nphi = length(phi);                 % number of $\phi$


