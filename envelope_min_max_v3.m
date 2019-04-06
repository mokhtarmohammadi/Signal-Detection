function [env_max, env_min, ex_no, zc_no, flag] = envelope_min_max_v3(a, ep, em); 
% Envelope calculator with 0 gradient spline interpolation
%
% This forms the initial stages of the EMD process
%
%
% Nathan Stevenson 
% SPR June 2004

N = length(a);
u = diff(a);
uu = a;

zc_max =[];
zc_min = [];

count_max = 1;
for zz = 1:(length(u)-1)
    if u(zz)>0
       if u(zz+1)<0
        zc_max(count_max) = zz;
        count_max = count_max+1;
    end
end
end

clear zz

count_min = 1;
for zz = 1:(length(u)-1)
    if u(zz)<0
       if u(zz+1)>0
        zc_min(count_min) = zz;
        count_min = count_min+1;
    end
end
end
clear zz

ex_no = length(zc_min) + length(zc_max);

zc_no = 0;
for zz = 1:(length(uu)-1)
    if uu(zz)>0
       if uu(zz+1)<0
        zc_no = zc_no+1;
    end
end
end

clear zz
for zz = 1:(length(uu)-1)
    if uu(zz)<0
       if uu(zz+1)>0
        zc_no = zc_no+1;
    end
end
end
clear zz

if length(zc_min)>=2 && length(zc_max)>=2;
% This section sets the end points of the cubic spline interpolation to the
% nearest maxima and minima value (an attempt to set the local mean of the 
% ends to zero)
        zc_max = zc_max+1;
        zc_min = zc_min+1;

        zcx = zc_max;
        acx = a(zc_max);
        ac_zcmax = a(zc_max);
        zzz1 = zc_max(2)+1;
        zc_max = [1 zzz1-zc_max(1) zc_max+zzz1];
        ac_zcmax = [ac_zcmax(2) ac_zcmax(1) ac_zcmax];
        zc_max = [zc_max 2*(N+zzz1)-zc_max(end)-1 2*(N+zzz1)-zc_max(end-1)-1];
        ac_zcmax = [ac_zcmax ac_zcmax(end) ac_zcmax(end-1)];
        zzz2 = max(zc_max)-zzz1-N;
        
        ac_zcmin = a(zc_min);
        zzz3 = zc_min(2)+1;
        zc_min = [1 zzz3-zc_min(1) zc_min+zzz3];
        ac_zcmin = [ac_zcmin(2) ac_zcmin(1) ac_zcmin];
        zc_min = [zc_min 2*(N+zzz3)-zc_min(end)-1 2*(N+zzz3)-zc_min(end-1)-1];
        ac_zcmin = [ac_zcmin ac_zcmin(end) ac_zcmin(end-1)];
        zzz4 = max(zc_min)-zzz3-N;

        
       z1 = 1:max(zc_max);
       z2 = 1:max(zc_min);
       env_max_1 = spline(zc_max, ac_zcmax, z1);
       env_min_1 = spline(zc_min, ac_zcmin, z2);         
       env_max = env_max_1(zzz1+1:end-zzz2);
       env_min = env_min_1(zzz3+1:end-zzz4);
       flag = 0;

else
    
env_min = ep;
env_max = em;
flag = 1;    

end


