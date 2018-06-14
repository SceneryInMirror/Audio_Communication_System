function [mpsk_sin, mpsk_cos] = mpsk(bits, Ns)

sequence_sin = [];
sequence_cos = [];
mpsk_sin = [];
mpsk_cos = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate Digital Signal %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
amp_sin = [1, 1/sqrt(2), 0, -1/sqrt(2), -1, -1/sqrt(2), 0, 1/sqrt(2)]; % amplitude options for Q_signal 
amp_cos = [0, 1/sqrt(2), 1, 1/sqrt(2), 0, -1/sqrt(2), -1, -1/sqrt(2)]; % amplitude options for I_signal
Nb = 100; % point number in digital signal per bit
dphi = 2 * pi / Nb;
phi = dphi:dphi:2*pi;
for n = 1:Ns
    k = 4 * bits(1, n) + 2 * bits(2, n) + bits(3, n) + 1;
    sequence_sin = [sequence_sin amp_sin(k) * ones(1, Nb)];
    sequence_cos = [sequence_cos amp_cos(k) * ones(1, Nb)];
end
    
figure(1)
subplot(2,2,1)
plot(sequence_sin);
subplot(2,2,2)
plot(sequence_cos);

%%%%%%%%%%%%%%%%%%%%
%%% Realize 8PSK %%%
%%%%%%%%%%%%%%%%%%%%
mpsk_sin = sequence_sin .* sin(repmat(phi, [1,Ns]));
mpsk_cos = sequence_cos .* cos(repmat(phi, [1,Ns]));
