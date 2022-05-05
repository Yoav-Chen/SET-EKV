
A=[];
for i=1:121
    A=[A,Ids(i,:)];
end

A=A'
A=abs(A)

A=A+1E-12;
