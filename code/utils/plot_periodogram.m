
function [] = plot_periodogram( Y ) 

T = length(Y); 
if mod(T,2) == 1
    T = T-1;
end
x = Y(1:T);

Fs = 1/0.72; % sampling frequency 1 kHz

t = 0:1/Fs:((T*0.72)-1/Fs); % time points to sample at 
 

figure(1); clf; 
plot(x); 


xdft = fft(x); % discrete fourier transform 
% get only half of the imaginary frequencies. first (zero frequncey) and 
% N/2 + 1th (Nyquist frequency) entry are unique 
% For N/2 - 2 frequencies, the FFT will have the complex conjugate frequncy
% as well 
xdft = xdft(1:T/2+1); 
psdx = (1/(Fs*T)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1); % multiply by two to account for positive and negative imaginary parts 
freq = 0:Fs/length(x):Fs/2; % define frequency domain 

figure(2); 
plot(freq,10*log10(psdx)); % decibal is 10*log10( power ratio ) 1dB is smallest perceptible change in sound to humans 
grid on
title('Periodogram Using FFT')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')