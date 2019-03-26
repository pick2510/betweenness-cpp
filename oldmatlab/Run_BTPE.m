clear all
i=0;
A=sprintf('RadiusSorted.txt');
R=load(A);
B=sprintf('ParticleID.txt');
ID=load(B);

for j=46800000%:5000:47000000;      
   
    i=i+1;
    ncon=0;
    
    C=sprintf('contact%d.txt',j);
    contact=load(C);    
    ncon = length(contact);
    npar = 7996;
    
    PE=zeros(7996,1);
    Between=zeros(7996,1);
    Connectivity= zeros (7996,7996); 
        
    %contact info    
    m=0; PEtot=0;   
    for m=1:1:ncon
            p1=0; p2=0; 
            rad1=0; rad2=0;
            Radstar=0; Ystar=0; Gstar=0; PEnor=0; PEtan=0; Ftan=0; Fnor=0; kn=0; kt=0;
            
            id=zeros(2,2);
            id(1,1)= contact(m,13);
            id(2,1)= contact(m,14);            
            p1 = find(R==id(1,1));
            p2 = find(R==id(2,1));
            rad1=R(p1,2);
            rad2=R(p2,2);
            
            %Connectivity Matrix
    
            Connectivity(p1,p2)=1;
            Connectivity(p2,p1)=1;
            
            %Potential Energy
            
            Radstar=(rad1*rad2)/(rad1+rad2);
            Ystar=(6.5e11)/(2*(1-.25^2));
            Gstar=(6.5e11)/(4*(2-.25)*(1+.25));
            
            kn=(4/3)*Ystar*sqrt(Radstar*contact(m,32));
            kt=8*Gstar*sqrt(Radstar*contact(m,32));
            
            Fnor= sqrt(contact(m,19)*contact(m,19)+contact(m,20)*contact(m,20)+contact(m,21)*contact(m,21));
            Ftan= sqrt(contact(m,22)*contact(m,22)+contact(m,23)*contact(m,23)+contact(m,24)*contact(m,24));            
            
            PEnor=0.5*(Fnor*Fnor/kn);
            PEtan=0.5*(Ftan*Ftan/kt);
            PEtot=PEnor+PEtan;
            
            PE(p1,1)=PE(p1,1)+0.5*PEtot;
            PE(p2,1)=PE(p2,1)+0.5*PEtot;
            
    end
    Between=transpose(betweenness_bin(Connectivity));    
    
    OutPE(1,i)=j;
    OutPE(2:7997,i)=PE;
    OutBetween(1,i)=j;
    OutBetween(2:7997,i)=Between;
    OutCoord(i,1)=j;
    OutCoord(i,2)=ncon;
end
dlmwrite('PE.txt', OutPE,'delimiter',' ');
dlmwrite('Between.txt', OutBetween,'delimiter',' ');
dlmwrite('Coord.txt', OutCoord,'delimiter',' ');
