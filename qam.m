clear, clc, close all
tic

    N=2048;                     %N=number of sample
    len=10;                    %repeat times
    C=[-3+3j,-3+1j,-3-3j,-3-1j,-1+3j,-1+1j,-1-3j,-1-1j,3+3j,3+1j,3-3j,3-1j,1+3j,1+1j,1-3j,1-1j];% C=16-QAM fomula of grey
    k=16;                       %k=scale of points
    % AWGN initialize
    Es=10;                      %symbol energy when in 16-qam and d=1
    M=16;                       %16QAM
    Eb=Es/log2(M);              %bit energy
    
for SNR=0:16;                   %SNR=0-16dB,step=1
    err=0;                      %Initialization err number

    for negentropie=1:len
%%
%sending

% Initialize random generator
RN=sum(100*clock);
RS=RandStream('mt19937ar','seed',RN);
RandStream.setGlobalStream(RS);

%order
bk=randi([0,15],1,N);
Yk=C(bk+1);

d=ifft(Yk);                     %IFFT
X=[d(N-k+1:N),d];               %CP
X0=sqrt(N)*X;                   %energy normalization

%%
%channel path : sigma for AWGN,cont from AWGN initialize
N0=Eb/(10^(SNR/10));
sigma=sqrt(N0/2);

%add noise,Y0=output
Y0=X0+sigma*(randn(1,N+k)+1j*randn(1,N+k));

%%
%%Receiver part1
Y1=Y0((k+1):(k+N));%remove CP
Y2=fft(Y1)/sqrt(N);%


% send Yy into 16-qam demodul
%16-qam demodulator
for i=1:2048;
for K=1:16
Dist(K)=abs(Y2(i)-C(K));
end
%%var(i)=find(min(Dist));
[Dis(i),var(i)]=min(Dist);
bk_new(i)=var(i)-1; %% bk_new is output of system
end

%%
%%Compare with signal to compute EBR
err_N=biterr(bk_new,bk);
%A=[1 2 1 2 1];B=[1 2 1 2 7];biterr(A,B)===answer=2
err=err+err_N;%Initialization err number
    end

%%
%%Receiver part2
   
    BER(SNR+1)=err/(len*N*4);
    str=sprintf('SNR=%0.1f,BER=%0.6f',SNR,BER(SNR+1));%printf BER and SNR in just one SNR after 1000times
    disp(str)
end
toc

%plot 16-qam BER in AWGN  
SNRdb=[0:16];
semilogy(SNRdb,BER,'ro-');
hold on;

%plot the theoretical 16-qam & 2-qam BER in AWGN
SNRlin=10.^(SNRdb/10);
Psm=2*(1-1/sqrt(M))*qfunc(sqrt(3*4*SNRlin/(M-1)));%M=16,so k=4
PM=1-(1-Psm).^2;
Pb16=PM/4; %k=log2(M)=4
Pb2=qfunc(sqrt(2*SNRlin));
semilogy(SNRdb,Pb16,SNRdb,Pb2);
grid on;
axis([0 40 10^-6 1]); %scale of axis
xlabel('SNR in dB');
ylabel('BER in dB');
title('BER in AWGN CHANNEL');
legend('simulative 16QAM BER','theoretical 16QAM BER','theoretical BPSK BER',4);