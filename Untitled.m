
clc
clear
close all
SampFreq = 128;
addpath('D:\tfsa_5-5\windows\win64_bin');

t = 0:1/SampFreq:1-1/SampFreq;



% COde for 2 component
code=2;

if code==2
 Sig2 = 1*exp(1i*(2*pi*(5*t))); %300t或者150t
Sig1 =1*exp(1i*(2*pi*(15*t +2*t.^3)));

Sig3 = exp(1i*(2*pi*(35*t +2*t.^3)));
Sig4 =1*exp(1i*(2*pi*(62*t -10*t.^3)));
SigA =1*Sig1  +Sig2+Sig3;%+Sig4;%+1*Sig2;
IF_O(:,1)=63*t.^2+0;

IF_O(:,2)=-45*t.^2+45;
IF_O(:,3)=-45*t.^2+62;
 
num=3;
elseif code==3
 Sig2 = 1*exp(1i*(2*pi*(30*t))); %300t或者150t
Sig1 =1*exp(1i*(2*pi*(0*t +21*t.^3)));
Sig1 =1*exp(1i*(2*pi*(0*t +10*t)));

Sig3 = exp(1i*(2*pi*(45*t -15*t.^3)));
Sig4 =1*exp(1i*(2*pi*(53*t -15*t.^3)));
SigA =1*Sig1  +Sig3+Sig4;%+1*Sig2;
IF_O(:,1)=63*t.^2+0;

IF_O(:,2)=-45*t.^2+45;
IF_O(:,3)=-45*t.^2+53;
 
num=3;
end
IF_O=IF_O/(SampFreq/2);

Fs=128;

[imf] = emd(real(SigA),3);
Spec=quadtfd(imf(1,:),length(imf)-1,1,'specx',45,'hamm');

for ii=2:3
Spec1=quadtfd(imf(ii,:),length(imf)-1,1,'specx',45,'hamm');
Spec=Spec+Spec1;

end
tfsapl(SigA,Spec)
Spec=quadtfd(SigA,length(imf)-1,1,'specx',45,'hamm');
figure;
tfsapl(SigA,Spec)
%z = hilbert(imf(1,:));
%instfreq = Fs/(2*pi)*diff(unwrap(angle(Sig2)));


%hs = hht(imf,128) ;
