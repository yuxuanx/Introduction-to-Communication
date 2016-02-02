function [ Ifinal, Qfinal ] = decisionMaking( mf_phase )
s = [(1 + 1i) (1 - 1i) (-1 + 1i) (-1 - 1i)]/sqrt(2);
Ifinal=zeros(1,length(mf_phase));
Qfinal=zeros(1,length(mf_phase));
for i=1:length(mf_phase)
    D1=norm(s(1)-mf_phase(i));%Calculate euclidean distance to each
    D2=norm(s(2)-mf_phase(i));%point of our constellation
    D3=norm(s(3)-mf_phase(i));
    D4=norm(s(4)-mf_phase(i));
    D=[D1 D2 D3 D4]; %Put all the distance in one vector
    [~, I]=min(D);   %Search for the index of the smallest value
    Ifinal(i)=sqrt(2)*real(s(I)); %And use this index to determine which symbol was send
    Qfinal(i)=sqrt(2)*imag(s(I));
end
end

