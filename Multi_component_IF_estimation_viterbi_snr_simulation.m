clc
clear
close all
SampFreq = 128;
addpath('D:\tfsa_5-5\windows\win64_bin');

t = 0:1/SampFreq:1-1/SampFreq;



% COde for 2 component
code=3;

if code==2
    Sig2 = 1*exp(1i*(2*pi*(30*t))); %300t或者150t
    Sig1 =1*exp(1i*(2*pi*(0*t +21*t.^3)));
    
    Sig3 = exp(1i*(2*pi*(45*t -15*t.^3)));
    Sig4 =1*exp(1i*(2*pi*(62*t -15*t.^3)));
    SigA =1*Sig1  +Sig3+Sig4;%+1*Sig2;
    IF_O(:,1)=63*t.^2+0;
    
    IF_O(:,2)=-45*t.^2+45;
    IF_O(:,3)=-45*t.^2+62;
    
    num=3;
elseif code==3
    Sig2 = 1*exp(1i*(2*pi*(30*t))); %300t或者150t
    Sig1 =1*exp(1i*(2*pi*(0*t +21*t.^3)));
    
    Sig3 = exp(1i*(2*pi*(45*t -15*t.^3)));
    Sig4 =1*exp(1i*(2*pi*(53*t -15*t.^3)));
    SigA =1*Sig1  +Sig3+Sig4;%+1*Sig2;
    IF_O(:,1)=63*t.^2+0;
    
    IF_O(:,2)=-45*t.^2+45;
    IF_O(:,3)=-45*t.^2+53;
    
    num=3;
end
IF_O=IF_O/(SampFreq/2);


% HADTFD BASED

iiii=0;
dis=1;
%Sig =1*Sig1+Sig3;


NS=2;
for kkk=0:1:5
    iiii=0;
    
    
    
    
    for snr=0:2:10
        iiii=iiii+1;
        
        for k1=1:NS
            Sig =SigA;
            
            Sig=awgn(Sig,snr,'measured');
            if kkk==0   %ADTFD+VITERBI
                [findex] = ADTFD_IF_estimation_viterbi_modified(Sig, num);
                
            elseif kkk==1 %SPEC+RPGP
                [fidexmult] = Proposed_IF_estimation_spec_RPGP(Sig, num,5);
                [findex,interset] = RPRG(fidexmult,15);
            elseif kkk==2 % ADTFD+RPGP
                [fidexmult] = Proposed_IF_estimation_shot_duration(Sig, num,5);
                [findex,interset] = RPRG(fidexmult,15);
            elseif kkk==3
                %SPEC+viterbi
                [findex] = Proposed_IF_estimation_spec_viterbi(Sig, num);
            elseif kkk==4 %ADTFD+old viterbi
                [findex] = ADTFD_IF_estimation_viterbi_old(Sig, num);
            else %EMD
                [findex] = IF_estimation_EMD(Sig, num,5);
            end
            msee=0.1*ones(1,num);
            IF=zeros(1,length(Sig));
            dis=0;
            
            
            for ii=1:num
                
                t=1:128;
                IF=findex(ii,:)/length(Sig);
                t=t(5:end-5);
                for i=1:num
                    c(i)=sum(abs(IF(t)'-IF_O(t,i)).^2);
                end
                [a1 b1]=min(c);
                if msee(b1)>=a1(1)/length(Sig)
                    msee(b1)=a1(1)/length(Sig);
                end
                if dis==1
                    figure;
                    plot(t,IF(t),'-',t,IF_O(t,b1),'d');
                end
            end
            msee1(k1)=mean(msee);
        end
        
        if kkk==0
            mse_adtfd_viterbi_modified(iiii)=mean(msee1);
        elseif kkk==1
            mse_Spec_rpgp(iiii)=mean(msee1);
        elseif kkk==2
            mse_adtfd_rpgp(iiii)=mean(msee1);
        elseif kkk==3
            mse_Spec_viterbi(iiii)=mean(msee1);
        elseif kkk==4
            mse_adtfd_viterbi_old(iiii)=mean(msee1);
        else
            mse_emd(iiii)=mean(msee1);
            
        end
    end
end


snr=0:2:10;
plot(snr,10*(log10(mse_adtfd_viterbi_modified(1:6))),'-.b+','linewidth',4);
hold on;
plot(snr,10*(log10(mse_Spec_rpgp(1:6))),'--md','linewidth',4);

hold on;
plot(snr, 10*(log10(mse_adtfd_rpgp(1:6))),'-rh','linewidth',4);

hold on;
plot(snr, 10*(log10(mse_Spec_viterbi(1:6))),'-.y+','linewidth',4);
hold on;
plot(snr, 10*(log10(mse_adtfd_viterbi_old(1:6))),'-.k+','linewidth',4);
hold on;
plot(snr, 10*(log10(mse_emd(1:6))),'-g+','linewidth',4);


%         spacing=4:13:17;
%
%         plot(snr,10*(log10(mse_adtfd_viterbi_modified(1:2))),'-.b+','linewidth',4);
%
%         hold on;
%         plot(snr,10*(log10(mse_Spec_rpgp(1:2))),'--md','linewidth',4);
%
%         hold on;
%         plot(snr, 10*(log10(mse_adtfd_rpgp(1:2))),'-rh','linewidth',4);
%
%         hold on;
%         plot(snr, 10*(log10(mse_Spec_viterbi(1:2))),'-.y+','linewidth',4);
%          hold on;
%         plot(snr, 10*(log10(mse_adtfd_viterbi_old(1:2))),'-.k+','linewidth',4);
%