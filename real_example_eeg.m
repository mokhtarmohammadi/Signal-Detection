clear all;
close all;
load seizure;
%addpath('D:\D\win64_bin\win64_bin');
addpath('D:\tfsa_5-5\windows\win64_bin');
%Sig=sez_dat(150,:);
% Sig=decimate(Sig,2);
% %Sig=diff(Sig);
% Sig=hilbert(Sig);
% 
% [findex] =ADTFD_IF_estimation_viterbi_modified_eeg(Sig(1:128), 3);
% figure; plot(findex')
% 
%  Sig=sez_dat(150,:);
%  Sig=decimate(Sig,2);
%  Sig=hilbert(Sig);
% % 
% t=0:1/16:8-1/16;
%  [findex] =ADTFD_IF_estimation_viterbi_modified_eeg(Sig(1:128), 3);
% figure; plot(t,8*findex'/128)
% 
%  axis([0 8 0 8]); 
%  xlabel('Time(s)');
%  ylabel('Instantaneous Frequency (Hz)');
 
 
 Sig=sez_dat(150,:);
 Sig=decimate(Sig,2);
 n=0:127;
 %Sig=Sig+0.035*cos(2*pi*(0.01*n+0.00001*n.^3));
  %Sig=Sig+0.1*cos(2*pi*(0.00*n+0.002*n.^2));
% %Sig=diff(Sig);
% 
 %[Spec,orienttfd]=HTFD_neww(Sig,2,15,54);
% imagesc(Spec)
t=0:1/16:8-1/16;
Sig=hilbert(Sig);
[findex] =ADTFD_IF_estimation_viterbi_modified_eeg(Sig(1:128), 3);
figure; plot(t,8*findex'/128)

 axis([0 8 0 8]); 
 xlabel('Time(s)');
 ylabel('Instantaneous Frequency (Hz)');
 
 %%% Results in noise%
 Sig=sez_dat(150,:);
 Sig=decimate(Sig,2);
 Sig=awgn(Sig,10,'measured');
 %Sig=Sig+0.035*cos(2*pi*(0.01*n+0.00001*n.^3));
  %Sig=Sig+0.1*cos(2*pi*(0.00*n+0.002*n.^2));
% %Sig=diff(Sig);
 Sig=hilbert(Sig);
% 
 %[Spec,orienttfd]=HTFD_neww(Sig,2,15,54);
% imagesc(Spec)
t=0:1/16:8-1/16;
 [findex] =ADTFD_IF_estimation_viterbi_modified_eeg(Sig(1:128), 3);
figure; plot(t,8*findex'/128)

 axis([0 8 0 8]); 
 xlabel('Time(s)');
 ylabel('Instantaneous Frequency (Hz)');
 
 
 %%% Results in case of intersecting components%
 Sig=sez_dat(150,:);
 Sig=decimate(Sig,2);
 n=0:127;
 Sig=Sig+0.035*cos(2*pi*(0.01*n+0.00001*n.^3));
  %Sig=Sig+0.1*cos(2*pi*(0.00*n+0.002*n.^2));
% %Sig=diff(Sig);
 Sig=hilbert(Sig);
% 
 %[Spec,orienttfd]=HTFD_neww(Sig,2,15,54);
% imagesc(Spec)
t=0:1/16:8-1/16;
 [findex] =ADTFD_IF_estimation_viterbi_modified_eeg(Sig(1:128), 4);
figure; plot(t,8*findex'/128)

 axis([0 8 0 8]); 
 xlabel('Time(s)');
 ylabel('Instantaneous Frequency (Hz)');
 
 

%load S;

%Sig=seiz_ds_05_20Hz(23,:);
%Sig=seiz_ds_05_20Hz(49,:);
% Sig=seiz_ds_05_20Hz(31,:);
% Sig=decimate(Sig,2);
% %Sig=diff(Sig);
% Sig=hilbert(Sig);
% % [Spec,orienttfd]=HTFD_neww(Sig,2,12,48);
% % tfsapl(Sig,Spec)
% 
% 
% [findex] =ADTFD_IF_estimation_viterbi_modified_eeg(Sig(1:128), 3);
% figure; plot(findex')


%  Sig=seiz_ds_05_20Hz(23,:);
% Sig=decimate(Sig,2);
% % %Sig=diff(Sig);
%  Sig=hilbert(Sig);
%  [Spec,orienttfd]=HTFD_neww(Sig,2,8,48);
%   tfsapl(Sig,Spec)
% % 
% % 
%  [findex] =ADTFD_IF_estimation_viterbi_modified_eeg(Sig(1:128), 4);
%  figure; plot(findex')