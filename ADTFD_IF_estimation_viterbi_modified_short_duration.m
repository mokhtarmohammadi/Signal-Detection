function [IA,IF] = ADTFD_IF_estimation_viterbi_modified_short_duration(Sig, Thresh,dis)
% Extract ridges for multi-component signals.
% In each iteration,the signal component associated with the extrated ridge is
% Sig is original Signal
% Threshold is ratio of energy vs earlier iteration
% Dis
if nargin==2
    dis=0;
elseif nargin==1
    Thresh=0.1;
    dis=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  output   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fidexmult: obtained ridge index for multi-component signals
Spec12=zeros(length(Sig),length(Sig));
[Spec,orienttfd]=HTFD_neww(Sig,2,20,64);%
if dis==1
    figure;tfsapl( real(Sig), Spec,'SampleFreq',16, 'GrayScale','on' );
end
Spec=Spec/max(Spec(:));
Spec(Spec<0.05)=0;
c = findridges_new_viterbi_adtfd(Spec,orienttfd);%(Spec,orienttfd,delta);

IFF=(c)/(2*length(Sig));

Phase=2*pi*filter(1,[1 -1],IFF);
s_dechirp=exp(-1i*Phase);

% For each sensor do the following steps

L=2;
%TF filtering for each sensor
s1 = Sig.*(s_dechirp);
s2=fftshift(fft(s1));
PPP=length(s2)/2;
s3=zeros(1,length(Sig));
s3(PPP-L:PPP+L)=s2(PPP-L:PPP+L);
s2(PPP-L:PPP+L)=0;
extr_Sig=ifft(ifftshift(s3)).*conj(s_dechirp);
s2=ifft(ifftshift(s2)).*conj(s_dechirp);

%Sig(iii)=Sig(iii)-extr_Sig(iii);
Sig=s2;%-extr_Sig(iii);
[Spec11,orienttfd]=HTFD_neww(extr_Sig,2,8,48);
Spec12=Spec12+Spec11;
prev_iter=mean(abs(s3.^2));
num=1;
IF1(num,:) = c;
    mean(abs(s3.^2))

while   mean(abs(s3.^2))>Thresh*max(prev_iter)  %Spec=tfr_stft_high(Sig);

    
    [Spec,orienttfd]=HTFD_neww(Sig,2,20,64);%
    %Spec=min(Spec1,Spec2);
    if dis==1
        figure;tfsapl( real(Sig), Spec,'SampleFreq',16, 'GrayScale','on' );
    end
    
    Spec=Spec/max(Spec(:));
    Spec(Spec<0.05)=0;
    c = findridges_new_viterbi_adtfd(Spec,orienttfd);%(Spec,orienttfd,delta);
    
    
    IFF=(c)/(2*length(Sig));
    
    Phase=2*pi*filter(1,[1 -1],IFF);
    s_dechirp=exp(-1i*Phase);
    
    % For each sensor do the following steps
    
    L=2;
    %TF filtering for each sensor
    s1 = Sig.*(s_dechirp);
    s2=fftshift(fft(s1));
    PPP=length(s2)/2;
    s3=zeros(1,length(Sig));
    s3(PPP-L:PPP+L)=s2(PPP-L:PPP+L);
    s2(PPP-L:PPP+L)=0;
    extr_Sig=ifft(ifftshift(s3)).*conj(s_dechirp);
    s2=ifft(ifftshift(s2)).*conj(s_dechirp);
    
    %Sig(iii)=Sig(iii)-extr_Sig(iii);
    Sig=s2;%-extr_Sig(iii);
    [Spec11,orienttfd]=HTFD_neww(extr_Sig,2,8,48);
    Spec12=Spec12+Spec11;
    num=num+1;
        IF1(num,:) = c;
    prev_iter(num)=mean(abs(s3.^2));
    
    % Add code for IA estimation
    %IA(i,:)=Spec(
        mean(abs(s3.^2))

end
if dis==1
    figure;tfsapl( real(Sig), Spec12,'SampleFreq',16, 'GrayScale','on' );
end
for i=1:num-1
    for ii=1:length(Sig)
        IF(i,ii)=IF1(i,ii);
        IA(i,ii)=Spec12(IF(i,ii),ii);
       
    end
end
end
