function [fidexmult] = IF_estimation_EMD(Sig, num,delta)
% IF estimation using EMD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  output   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fidexmult: obtained ridge index for multi-component signals
[imf] = emd(real(Sig),3);

for i = 1:num
Spec=quadtfd(imf(i,:),length(imf(i,:))-1,1,'specx',45,'hamm');
%Spec=quadtfd(Sig,length(Sig)/4-1,1,'mb',0.05,128);

c = findridges(Spec,delta);%(Spec,orienttfd,delta);
%c = findridges_neww(Spec,orienttfd,delta);


 %IF=(c)/(2*length(Sig));

             
fidexmult(i,:) = c;

end

end