function IF = findridges_new_viterbi(Spec)
Path=zeros(size(Spec));
Pathweight=zeros(1,length(Spec));
Path(:,1)=1:length(Path);
for i=2:length(Spec)
    for ii=1:length(Spec) % all members of
        Pathw=zeros(1,length(Spec));
        % For transition from ii to j.
        
        for jj=1:length(Spec)
            rr=abs(jj-ii)-2;
            if rr<0
                rr=0;
                
            end
            
%             if angg<5
%                 angg=0;
%             end
             
            if i>2
            angg=Path(jj,i-1)-jj-(jj-ii);
            if (Path(jj,i-1)-jj)*(jj-ii)>0
                angg=0;
            else
                angg=8*abs(angg);
            end
            else
                angg=0;
            end
            Pathw(jj)=rank(Spec(:,i),ii)+Pathweight(jj)+8*rr+angg;% weight of ii +
                   %     Pathw(jj)=rank(Spec(:,i),ii)+Pathweight(jj)+5*rr+angg;% weight of ii +

        end
        [value,index]=min(Pathw);
        Pathweight1(ii)=value;
        Path(ii,i)=index;
    end
    Pathweight=Pathweight1;
    
end
IF(1,length(Spec))=0;
[~,index]=min(Pathweight);
for i=length(Spec):-1:1
    IF(i)=index;
    index=Path(index,i);
end
% Code for decoding
end
function a= rank(A,b)
As=sort(A,'descend');
ind=find(As==A(b));
if A(b)~=0
a=ind(1)-1;
else
    a=128;

end
end