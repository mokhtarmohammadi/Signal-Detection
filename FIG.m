close all;
clear all;
fs=200;
M=1;
n=0:127;
N_S=10000;
addpath('D:\tfsa_5-5\windows\win64_bin');
P=0;
Q=0;
llll=0;
Ruho          = 2.0; % Noise Power Uncertainity Parameter 1 no ucertainity and 2 maximum power uncertainity.

s=exp(2*pi*1i*(0.48*n-0*0.2*n.^2/(2*128)-1.5*0.4*n.^3/(128*128*3)));%+exp(2*pi*1i*(0.3*n-0*0.2*n.^2/(2*128)-0.3*n.^3/(128*128*3)));
%s=exp(2*pi*1i*(0.45*n-0*0.2*n.^2/(2*128)-0.5*n.^3/(128*128*3)));%+exp(2*pi*1i*(0.3*n-0*0.2*n.^2/(2*128)-0.3*n.^3/(128*128*3)));
%s=exp(2*pi*1i*(0.5*n-0.5*n.^2/(2*128)));%+0.5*exp(2*pi*1i*(0.25*n-0.2*n.^2/(2*128)));%+exp(2*pi*1i*(0.3*n-0*0.2*n.^2/(2*128)-0.3*n.^3/(128*128*3)));
T1 = zeros(1,N_S);
T2 = zeros(1,N_S);
T3 = zeros(1,N_S);
T4 = zeros(1,N_S);
for SNR=-10:2:4
    llll=llll+1;
    
    s=exp(2*pi*1i*(0.48*n-0*0.2*n.^2/(2*128)-1.5*0.4*n.^3/(128*128*3)));%+exp(2*pi*1i*(0.3*n-0*0.2*n.^2/(2*128)-0.3*n.^3/(128*128*3)));
    P=0;
    Q=0;
    for i=1:N_S
        sorig = s;
        
        X = s;                             % mixed source
        ii = sign(randn);
        if  ii > 0
            Mul = 1;
            c(i) = 1;
            P = P + 1;
        else
            Mul = 0;
            c(i) = 0;
            Q = Q + 1;
        end
        %----------------------------------------------------
        %----------------------------------------------------
        sigma = 10^(-SNR/20);
        StdN = sqrt(sigma);
        % generate noise
        
        NoisePower_Uncetainity   = 1/Ruho + (Ruho - (1/Ruho)).*rand;
        Noise_UNC = NoisePower_Uncetainity*ones(M,1);
        StdN_NU   = StdN.*NoisePower_Uncetainity;
        
        
        w = StdN_NU.*(randn(M,128) + 1j*(randn(M,128)))./sqrt(2); % noise
        X=Mul*X+w;
        %----------------------------------------------------
        %----------------------------------------------------
        %X(randsample(128,32))=0;
        %X(2:3:end)=0;
        
        %X=w;
        %Summation of Auto Wigner
       % I=quadtfd(X(1,:),length(s)/2-1,1,'mb',0.2,length(X));%+quadtfd(X(2,:),length(s)/2-1,1,'mb',0.05,length(X))+quadtfd(X(3,:),length(s)/2-1,1,'mb',0.05,length(X));%+quadtfd(X(4,:),length(s)/2-1,1,'mb',0.05,length(X));
        %I=HTFD_neww(X,3,20,64);
        %I=quadtfd(s,length(s)/4-1,1,'mb',0.05,length(X));%+quadtfd(X(2,:),length(s)/2-1,1,'mb',0.05,length(X))+quadtfd(X(3,:),length(s)/2-1,1,'mb',0.05,length(X));%+quadtfd(X(4,:),length(s)/2-1,1,'mb',0.05,length(X));
        %[IA,IF]=ADTFD_IF_estimation_viterbi_modified_short_duration(X, 0.99,0);
        [IA,IF]=ADTFD_IF_estimation_viterbi_modified_IA(X, 1);
        %[IA,IF]=max(I);
        %IF=IF_COMPUTE(I);
        IF=IF(1,:);
        IF=IF(:);
        IA=IA(1,:);
        IA=IA(:);
        
        IF=IF-1;
        IF=IF/(2*length(s));
        
        %l=zeros(size(I));
        
        
        Phase=2*pi*filter(1,[1 -1],IF);
        sig_den=X.*exp(-1i*Phase.');
        %  sig_den=s.*exp(-1i*Phase);
        %I=HTFD_new3(X,3,20,64);
        NN=5;
        %NN=7;
        J=xcorr(sig_den,NN);
        
        
        
        %NN=7;
        J=xcorr(sig_den,NN);
        
        AA=toeplitz(J(NN+1:end));
        
        T1(i)=abs(prod(diag(AA)))/abs(det(AA));
        %T1(i)=sum(abs(AA(:)))/sum(abs(diag(AA)));
        %T1(i)=-abs(J(NN+1));
        
        %T3(i)=-sum(IA);
        %T3(i)=-sum(IA);
        T4(i)=sum(abs(X.^2));
        T3(i)=sum(abs(IA));
        %AA=toeplitz(xcorr(X,5));
        NN=5;
        J=xcorr((X),NN);
        AA=toeplitz(J(NN+1:end));
        %T2(i) =abs(det(AA))/abs(prod(diag(AA)));
        T2(i) =abs(prod(diag(AA)))/abs(det(AA));
        
        
        I=quadtfd(X(1,:),length(s)/2-1,1,'mb',0.2,length(X));
        [IA,IF]=max(I);
        IF=IF-1;
        IF=IF/(2*length(s));
        
        Phase=2*pi*filter(1,[1 -1],IF);
        sig_den=X.*exp(-1i*Phase);
        NN=5;
        J=xcorr(sig_den,NN);
        
        AA=toeplitz(J(NN+1:end));
        T6(i)=abs(prod(diag(AA)))/abs(det(AA));
        
        
        %T2(i) =sum(abs(AA(:)))/sum(abs(diag(AA)));
        %T2(i)=-sum(abs(X.^2));
        %T2(i)=sum(abs(AA))/sum(abs(diag(AA)));
        
    end
    [PF1, PD1] = roc1(T1,c,Q,P);
    %trapz(PD1,PF1)
    
    [PF2, PD2] = roc1(T2,c,Q,P);
    %trapz(PD2,PF2)
    
    [PF3, PD3] = roc1(T3,c,Q,P);
    
    [PF4, PD4] = roc1(T4,c,Q,P);
    
    [PF5, PD5] = roc1(T6,c,Q,P);
    
    %trapz(PD3,PF3)
    
    for iii=1:length(PF1)
        if PF1(iii+1)>=0.01
            PDD1(llll)= PD1(iii);
            break;
        end
    end
    for iii=1:length(PF2)
        if PF2(iii+1)>=0.01
            PDD2(llll)= PD2(iii);
            break;
        end
    end
    for iii=1:length(PF3)
        if PF3(iii+1)>=0.01
            PDD3(llll)= PD3(iii);
            break;
        end
    end
    for iii=1:length(PF4)
        if PF4(iii+1)>=0.01
            PDD4(llll)= PD4(iii);
            break;
        end
    end
    
    for iii=1:length(PF5)
        if PF5(iii+1)>=0.01
            PDD5(llll)= PD5(iii);
            break;
        end
    end
end
SNR=-10:2:4;
%SNR=-6:2:-4;
plot(SNR,PDD1,'o-');
hold on;
plot(SNR,PDD2,'rx-');
hold on;
plot(SNR,PDD3,'y*-');
hold on;
plot(SNR,PDD4,'k--');
hold on;
plot(SNR,PDD5,'g--');

xlabel('Signal to noise ratio (dB)');
ylabel('Detection Probability');
legend('The Proposed method', 'Covariance based detection','Time-frequency ridge energy','Energy detection','Conventional IF method');
PDD1
PDD2
PDD3
%legend('The Proposed method', 'time-domain covariance detector','Energy detection');