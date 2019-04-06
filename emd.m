function [imf, res] = emd(y, emd_no);
% Empirical Mode Decomposition - this function performs an empirical mode
% decomposition on the input signal.
%
% INPUT:- 
%         y - signal under analysis of length N
%         emd_no - specifies the specific function that fits the signal
%                  envelope (numeric 3 or 5 at this stage)
%
% OUTPUTS:-
%         imf - a matrix contain the intrinsic mode functions of the EMD
%               this matix is LxN where L is the number of IMFs
%         res - the residual
%
% Functions called
%    envelope_min_max -> fits the maxima and minima signal envelopes using
%                        cubic splines
%
% Nathan Stevenson
% SPR June 2005

x = y;
imf_no = 0;
zc_no = 5;
ex_no = 6;
env_max = 1; env_min = 1; flag = 0;
while (zc_no > 1) && (ex_no > 1) && (flag == 0)

    imf_no = imf_no+1;
    
if imf_no == 1
    yy = y; 
else
    yy = res;
end
  
clear mean_env
yyy = yy;    
count = 0;
while (count<1) && (flag == 0)
    eval(['[env_max, env_min, ex_no, zc_no, flag] = envelope_min_max_v' num2str(emd_no) '(yy, env_max, env_min);']);
    mean_env = (env_max+env_min)/2;
    yy = yy-mean_env;
%    ex_no-zc_no
    if ex_no==(zc_no-1) || (ex_no==zc_no+1) || (ex_no==zc_no)
        count = count+1;
    else
        count=0;
    end
end

imf(imf_no,:) = yy;
res = yyy-yy;
clear yy yyy

end