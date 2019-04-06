clear all;
close all;
load seizure;
%addpath('D:\D\win64_bin\win64_bin');
SampFreq = 128;
t = 0:1/SampFreq:1-1/SampFreq;
IF_O(:,1)=45*t.^2+10;
    
    IF_O(:,2)=-45*t.^2+45;
    IF_O(:,3)=-45*t.^2+62;
    IF_O=IF_O/(SampFreq/2);

addpath('D:\tfsa_5-5\windows\win64_bin');
Sig2 = 1*exp(1i*(2*pi*(30*t))); %300tªÚ’ﬂ150t
    Sig1 =1*exp(1i*(2*pi*(10*t +15*t.^3)));
    
    Sig3 = exp(1i*(2*pi*(45*t -15*t.^3)));
    %Sig3=[ zeros(1,32) Sig3(33:128-32) zeros(1,32)];
    Sig4 =1*exp(1i*(2*pi*(55*t -15*t.^3)));
    %Sig1=[ zeros(1,32) Sig1(33:128-32) zeros(1,32)];
    %Sig3=[ zeros(1,32) Sig3(33:128-32) zeros(1,32)];

    Sig =2*Sig1  +Sig3+0.5*Sig4;%+1*Sig2;
    Sig=awgn(Sig,1);
[IA,findex] =ADTFD_IF_estimation_viterbi_modified_short_duration(Sig(1:128), 0.05,1);
findex=findex/(SampFreq);
  figure; plot(findex')
  hold on;plot(IF_O,':');