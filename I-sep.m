%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Field Dimensions - x and y maximum (in meters)
xm=100;
ym=100;
%x and y Coordinates of the Sink
sink.x=0.5*xm;
sink.y=0.5*ym;
%Number of Nodes in the field
n=100
allive1=n;
allive2=n;
allive3=n;
allive7=n;
%Optimal Election Probability of a node
%to become cluster head
p=0.1;
sv=0;                                  %%%%%%previously Sensed value S(v)
b=0.5;
ph=0.5;
%Energy Model (all values in Joules)
%Initial Energy 
Eo=0.5;
%Eelec=Etx=Erx
ETX=50*0.000000001;
ERX=50*0.000000001;
%Transmit Amplifier types
Efs=10*0.000000000001;
Emp=0.0013*0.000000000001;
%Data Aggregation Energy
EDA=5*0.000000001;
h=100;                                 %%%%%%Hard Thres%%%%%hold H(t)
s=2;                                   %%%%%%Soft thres%%%%%hold  S(t)
%temprature range
tempi=50;
tempf=200;
%Values for Hetereogeneity
%Percentage of nodes than are advanced
m=0.1;
%\alpha
a=1;

u=a/2;
%maximum number of rounds
rmax=8000
allive=n;
E_adv=Eo*(1+a);
%%%%%%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%
packets_TO_BS=0;
%Computation of do
do=sqrt(Efs/Emp); 
rcountCHs1=0;
rcountCHs2=0;
rcountCHs3=0;
rcountCHs4=0;
rcountCHs5=0;
rcountCHs6=0;
rcountCHs7=0;

for i=1:1:n
    S1(i).xd=rand(1,1)*xm;
    S2(i).xd=S1(i).xd;
    S3(i).xd=S1(i).xd;
    S4(i).xd=S1(i).xd;
    S6(i).xd=S1(i).xd;
    S7(i).xd=S1(i).xd;
    XR7(i)=S7(i).xd;
    
    XR6(i)=S6(i).xd;
    
    XR4(i)=S4(i).xd;
    XR3(i)=S3(i).xd;
    XR2(i)=S2(i).xd;
    XR1(i)=S1(i).xd;
    S1(i).yd=rand(1,1)*ym;
    S2(i).yd=S1(i).yd;
    S3(i).yd=S1(i).yd;
    S4(i).yd=S1(i).yd;
    S6(i).yd=S1(i).yd;
    S7(i).yd=S1(i).yd;
     YR7(i)=S7(i).yd;
     YR6(i)=S6(i).yd;
     S7(i).G=0;
    S6(i).G=0;
    YR4(i)=S4(i).yd;
    S4(i).G=0;
    YR3(i)=S3(i).yd;
    S3(i).G=0;
    YR2(i)=S2(i).yd;
    YR1(i)=S1(i).yd;
    S1(i).G=0;
    S2(i).G=0;
    
    %%%%%%%

    %initially there are no cluster heads only nodes
    S1(i).type='N';
    S2(i).type='N';
    S3(i).type='N';
    S4(i).type='N';
     S6(i).type='N';
     S7(i).type='N';
     S6(i).checked=0;
     
    temp_rnd0=i;
   
    %Random Election of Normal Nodes
    if (temp_rnd0>=m*n+1) 
        S1(i).E=Eo;
        S1(i).ENERGY=0;
        %%%%%plot(S(i).xd,S(i).yd,'o');
        hold on;
    end
    %Random Election of Advanced Nodes
    if (temp_rnd0<m*n+1)  
        S1(i).E=Eo*(1+a)
        S1(i).ENERGY=1;
        %%%%%%plot(S(i).xd,S(i).yd,'+');
         hold on;
    end
    
     temp_rnd0=i;
    %Random Election of Normal Nodes
    if (temp_rnd0>=m*n+1) 
        S2(i).E=Eo;
        S2(i).ENERGY=0;
        %%%%%%plot(S(i).xd,S(i).yd,'o');
        hold on;
    end
    %Random Election of Advanced Nodes
    if (temp_rnd0<m*n+1)  
        S2(i).E=Eo*(1+a)
        S2(i).ENERGY=1;
        %%%%%%plot(S(i).xd,S(i).yd,'+');
         hold on;
    end
   
    figure(1);
    temp_rnd0=i;
    %Random Election of Normal Nodes
    if (temp_rnd0>=m*n+1&& temp_rnd0<50) 
        S3(i).E=Eo;
        S3(i).ENERGY=0;
        plot(S3(i).xd,S3(i).yd,'r.','LineWidth',2);
        %hold on;
    end
    %Random Election of Advanced Nodes
    if (temp_rnd0<m*n+1)  
        S3(i).E=Eo*(1+a)
        S3(i).ENERGY=1;
        plot(S3(i).xd,S3(i).yd,'b+','LineWidth',2);
         %hold on;
    end
    %Random Election of intermediate Nodes
    if (temp_rnd0>=50)  
        S3(i).E=Eo*(1+u)
        S3(i).ENERGY=1.5;
        plot(S3(i).xd,S3(i).yd,'g*','LineWidth',2);
    end
   
    temp_rnd0=i;
    %Random Election of Normal Nodes
    if (temp_rnd0>=m*n+1) 
        S4(i).E=Eo;
        S4(i).ENERGY=0;
        %%%%%%plot(S(i).xd,S(i).yd,'o');
        hold on;
    end
    %Random Election of Advanced Nodes
    if (temp_rnd0<m*n+1)  
        S4(i).E=Eo*(1+a)
        S4(i).ENERGY=1;
        %%%%%%plot(S(i).xd,S(i).yd,'+');
         hold on;
    end
    temp_rnd0=i;
    %Random Election of Normal Nodes
    if (temp_rnd0>=m*n+1) 
        S6(i).E=Eo;
        S6(i).ENERGY=0;
        %%%%%%plot(S(i).xd,S(i).yd,'o');
        hold on;
    end
    %Random Election of Advanced Nodes
    if (temp_rnd0<m*n+1)  
        S6(i).E=Eo*(1+a)
        S6(i).ENERGY=1;
        %%%%%%plot(S(i).xd,S(i).yd,'+');
         hold on;
    end
    
    temp_rnd0=i;
    %Random Election of Normal Nodes
    if (temp_rnd0>=m*n+1 && temp_rnd0<50) 
        S7(i).E=Eo;
        S7(i).ENERGY=0;
        %plot(S7(i).xd,S7(i).yd,'o');
        hold on;
    end
    %Random Election of Advanced Nodes
    if (temp_rnd0<m*n+1)  
        S7(i).E=Eo*(1+a)
        S7(i).ENERGY=1;
        %plot(S7(i).xd,S7(i).yd,'+');
         hold on;
    end
     %Random Election of intermediate Nodes
    if (temp_rnd0>=50)  
        S7(i).E=Eo*(1+u)
        S7(i).ENERGY=1.5;
        %plot(S7(i).xd,S7(i).yd,'-');
    end
   
end



S1(n+1).xd=sink.x; %sink is a n+1 node, x-axis postion of a node
S1(n+1).yd=sink.y; %sink is a n+1 node, y-axis postion of a node
S2(n+1).xd=sink.x; %sink is a n+1 node, x-axis postion of a node
S2(n+1).yd=sink.y; %sink is a n+1 node, y-axis postion of a node
S3(n+1).xd=sink.x; %sink is a n+1 node, x-axis postion of a node
S3(n+1).yd=sink.y; %sink is a n+1 node, y-axis postion of a node
S4(n+1).xd=sink.x; %sink is a n+1 node, x-axis postion of a node
S4(n+1).yd=sink.y; %sink is a n+1 node, y-axis postion of a node
S6(n+1).xd=sink.x; %sink is a n+1 node, x-axis postion of a node
S6(n+1).yd=sink.y; %sink is a n+1 node, y-axis postion of a node

S7(n+1).xd=sink.x;
S7(n+1).yd=sink.y;
countCHs1=0;  %variable, counts the cluster head
countCHs2=0;  %variable, counts the cluster head
countCHs3=0;  %variable, counts the cluster head
countCHs4=0;  %variable, counts the cluster head
countCHs6=0;  %variable, counts the cluster head
countCHs7=0;
%counter for CHs per round
rcountCHs7=0;
cluster7=1;
cluster6=1;  %cluster is initialized as 1
cluster1=1;  %cluster is initialized as 1
cluster2=1;  %cluster is initialized as 1
cluster3=1;  %cluster is initialized as 1
cluster4=1;  %cluster is initialized as 1
flag_first_dead1=0; %flag tells the first node dead
flag_first_dead2=0; %flag tells the first node dead
flag_first_dead3=0; %flag tells the first node dead
flag_first_dead4=0; %flag tells the first node dead
flag_first_dead6=0; %flag tells the first node dead
flag_first_dead7=0;


dead1=0;  %dead nodes count initialized to 0
dead2=0;  %dead nodes count initialized to 0
dead3=0;  %dead nodes count initialized to 0
dead4=0;  %dead nodes count initialized to 0
dead6=0;  %dead nodes count initialized to 0
dead7=0;  %dead nodes count initialized to 0
first_dead6=0;
first_dead1=0;
first_dead2=0;
first_dead3=0;
first_dead4=0;
first_dead7=0;
allive1=n;
allive2=n;
allive3=n;
allive4=n;
allive6=n;
allive7=n;
%counter for bit transmitted to Bases Station and to Cluster Heads
packets_TO_BS1=0;
packets_TO_BS2=0;
packets_TO_BS3=0;
packets_TO_BS4=0;
packets_TO_BS6=0;
packets_TO_BS7=0;
packets_TO_CH7=0;
packets_TO_CH6=0;
packets_TO_CH1=0;
packets_TO_CH2=0;
packets_TO_CH3=0;
packets_TO_CH4=0;

for r=0:1:rmax
  
    r
 cv = tempi + (tempf-tempi).*rand(1,1);  %%%%%%Current sensing value C(v)
  %Election Probability for Normal Nodes
  pnrm3=( p/ (1+a*m+b*u) );
  %Election Probability for Advanced Nodes
  padv3= ( p*(1+a)/(1+a*m+b*u) );
   %Election Probability for intermediate Nodes
  pint3= ( p*(1+u)/(1+a*m+b*u) );
    %Operation for heterogeneous epoch
  if(mod(r, round(1/pnrm3) )==0)
    for i=1:1:n
        S3(i).G=0;
        S3(i).cl=0;
    end
  end
  %Operation for heterogeneous intermediate nodes epoch
  if(mod(r, round(1/pint3) )==0)
    for i=1:1:n
        if (S3(i).ENERGY==1.5)
        S3(i).G=0;
        S3(i).cl=0;
    end
    end
  end

 %Operations for advance nodes sub-epochs
 if(mod(r, round(1/padv3) )==0)
    for i=1:1:n
        if(S3(i).ENERGY==1)
            S3(i).G=0;
            S3(i).cl=0;
        end
    end
  end

 hold off;

%Number of dead nodes
dead3=0;
%Number of dead Advanced Nodes
dead_a3=0;
%Number of dead Normal Nodes
dead_n3=0;
%Number of dead intermediate Nodes
dead_I3=0;

%counter for bit transmitted to Bases Station and to Cluster Heads

packets_TO_CH3=0;
%counter for bit transmitted to Bases Station and to Cluster Heads 
%per round
PACKETS_TO_CH3(r+1)=0;


figure(1);

for i=1:1:n
    %checking if there is a dead node
    if (S3(i).E<=0)
        plot(S3(i).xd,S3(i).yd,'k.','LineWidth',2);
        title('Dead nodes');
        dead3=dead3+1;
        if(S3(i).ENERGY==1)
            dead_a3=dead_a3+1;
        end
        if(S3(i).ENERGY==1.5)
            dead_I3=dead_I3+1;
        end
        if(S3(i).ENERGY==0)
            dead_n3=dead_n3+1;
        end
        hold on;    
    end
    if S3(i).E>0
        S3(i).type='N';
        if (S3(i).ENERGY==0)  
        plot(S3(i).xd,S3(i).yd,'bo','LineWidth',2);
        end
        if (S3(i).ENERGY==1.5)  
        plot(S3(i).xd,S3(i).yd,'r*','LineWidth',2);
        end
        if (S3(i).ENERGY==1)  
        plot(S3(i).xd,S3(i).yd,'g+','LineWidth',2);
        end
        hold on;
    end
end
%plot(S(n+1).xd,S(n+1).yd,'x');


STATISTICS(r+1).DEAD3=dead3;
DEAD3(r+1)=dead3;
%ALLIVE(r+1)=(allive-dead);
DEAD_N3(r+1)=dead_n3;
DEAD_A3(r+1)=dead_a3;
DEAD_I3(r+1)=dead_I3;

ALIVE3(r+1)=allive3-dead3;
%When the first node dies
if (dead3==1)
    if(flag_first_dead3==0)
        first_dead3=r
        flag_first_dead3=1;
    end
end

countCHs3=0;
cluster3=1;
for i=1:1:n
   if(S3(i).E>0)
   temp_rand3=rand;     
   if ( (S3(i).G)<=0)

 %Election of Cluster Heads for normal nodes
 if( ( S3(i).ENERGY==0 && ( temp_rand3 <= ( pnrm3 / ( 1 - pnrm3 * mod(r,round(1/pnrm3)) )) ) )  )

            countCHs3=countCHs3+1;
            
            packets_TO_BS3=packets_TO_BS3+1;
            PACKETS_TO_BS3(r+1)=packets_TO_BS3;
            
            S3(i).type='C';
            S3(i).G=100;
            C3(cluster3).xd=S3(i).xd;
            C3(cluster3).yd=S3(i).yd;
            plot(S3(i).xd,S3(i).yd,'k>','LineWidth',2);
            
            distance3=sqrt( (S3(i).xd-(S3(n+1).xd) )^2 + (S3(i).yd-(S3(n+1).yd) )^2 );
            C3(cluster3).distance3=distance3;
            C3(cluster3).id=i;
            X(cluster3)=S3(i).xd;
            Y(cluster3)=S3(i).yd;
            cluster3=cluster3+1;
            
            %Calculation of Energy dissipated
              distance3;
          if (cv >= h)
            test = cv-sv;
            if (test >= s)
            if (distance3>do)
                S3(i).E=S3(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance3*distance3*distance3*distance3 )); 
            end
            if (distance3<=do)
                S3(i).E=S3(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance3 * distance3 )); 
            end
            
            end
          end
            
 end
 %Election of Cluster Heads for intermediate  nodes
 if( ( S3(i).ENERGY==1.5 && ( temp_rand3 <= ( pint3 / ( 1 - pint3 * mod(r,round(1/pint3)) )) ) )  )

            countCHs3=countCHs3+1;
            packets_TO_BS3=packets_TO_BS3+1;
            PACKETS_TO_BS3(r+1)=packets_TO_BS3;
            
            S3(i).type='C';
            S3(i).G=100;
            C3(cluster3).xd=S3(i).xd;
            C3(cluster3).yd=S3(i).yd;
            plot(S3(i).xd,S3(i).yd,'k>','LineWidth',2);
            
            distance3=sqrt( (S3(i).xd-(S3(n+1).xd) )^2 + (S3(i).yd-(S3(n+1).yd) )^2 );
            C3(cluster3).distance3=distance3;
            C3(cluster3).id=i;
            X(cluster3)=S3(i).xd;
            Y(cluster3)=S3(i).yd;
            cluster3=cluster3+1;
            
            %Calculation of Energy dissipated
              distance3;
          if (cv >= h)
            test = cv-sv;
            if (test >= s)
                distance3;
            if (distance3>do)
                S3(i).E=S3(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance3*distance3*distance3*distance3 )); 
            end
            if (distance3<=do)
                S3(i).E=S3(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance3 * distance3 )); 
            end
            end     
          end
 end
 
  
        
 %Election of Cluster Heads for Advanced nodes
 if( ( S3(i).ENERGY==1 && ( temp_rand3 <= ( padv3 / ( 1 - padv3 * mod(r,round(1/padv3)) )) ) )  )
        
            countCHs3=countCHs3+1;
            packets_TO_BS3=packets_TO_BS3+1;
            PACKETS_TO_BS3(r+1)=packets_TO_BS3;
            
            S3(i).type='C';
            S3(i).G=100;
            C3(cluster3).xd=S3(i).xd;
            C3(cluster3).yd=S3(i).yd;
            plot(S3(i).xd,S3(i).yd,'k>','LineWidth',2);
            
            distance3=sqrt( (S3(i).xd-(S3(n+1).xd) )^2 + (S3(i).yd-(S3(n+1).yd) )^2 );
            C3(cluster3).distance3=distance3;
            C3(cluster3).id=i;
            X(cluster3)=S3(i).xd;
            Y(cluster3)=S3(i).yd;
            cluster3=cluster3+1;
            
            %Calculation of Energy dissipated
              distance3;
          if (cv >= h)
            test = cv-sv;
            if (test >= s)
            if (distance3>do)
                S3(i).E=S3(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance3*distance3*distance3*distance3 )); 
            end
            if (distance3<=do)
                S3(i).E=S3(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance3 * distance3 )); 
            end
            
            end
          end
        end     
    
    end
  end 
end

STATISTICS3(r+1).CLUSTERHEADS3=cluster3-1;
CLUSTERHS3(r+1)=cluster3-1;



%Election of Associated Cluster Head for Normal Nodes
for i=1:1:n
   if ( S3(i).type=='N' && S3(i).E>0 )
     if(cluster3-1>=1)
       min_dis3=sqrt( (S3(i).xd-S3(n+1).xd)^2 + (S3(i).yd-S3(n+1).yd)^2 );
       min_dis_cluster3=1;
       for c3=1:1:cluster3-1
           temp3=min(min_dis3,sqrt( (S3(i).xd-C3(c3).xd)^2 + (S3(i).yd-C3(c3).yd)^2 ) );
           if ( temp3<min_dis3 )
               min_dis3=temp3;
               min_dis_cluster3=c3;
           end
       end
       
       %Energy dissipated by associated Cluster Head
            min_dis3;
              
          if (cv >= h)
            test = cv-sv;
            if (test >= s)
            if (min_dis3>do)
                S3(i).E=S3(i).E- ( ETX*(4000) + Emp*4000*( min_dis3 * min_dis3 * min_dis3 * min_dis3)); 
            end
            if (min_dis3<=do)
                S3(i).E=S3(i).E- ( ETX*(4000) + Efs*4000*( min_dis3 * min_dis3)); 
            end
        %Energy dissipated
        if(min_dis3>0)
            S3(C3(min_dis_cluster3).id).E = S3(C3(min_dis_cluster3).id).E- ( (ERX + EDA)*4000 ); 
         PACKETS_TO_CH3(r+1)=n-dead3-cluster3+1; 
        end
            end
          end
       S3(i).min_dis3=min_dis3;
       S3(i).min_dis_cluster3=min_dis_cluster3;
           
   end
 end
end
%hold on;

countCHs3;
rcountCHs3=rcountCHs3+countCHs3;
sv=cv;
inisv=sv;
if (cv >= 180 && r<=2500)
            sv=inisv;
      end 

%Code for Voronoi Cells
%Unfortynately if there is a small
%number of cells, Matlab's voronoi
%procedure has some problems

%[vx,vy]=voronoi(X,Y);
%%plot(X,Y,'r*',vx,vy,'b-');
% %hold on;
% voronoi(X,Y);
% axis([0 xm 0 ym]);
%STATISTICS.ALLIVE(r+1)
P3.PACKETS_TO_BS3(r+1)=packets_TO_BS3;

end
r=0:8000;
figure(2);
plot(r,DEAD3,'--b');
legend('I-SEP');
xlabel('Number of rounds');
ylabel('Dead nodes');
title('Nodes dead during rounds');

figure(3);
%subplot(2,2,2);
plot(r,ALIVE3,'--r');
legend('I-SEP');
xlabel('Number of rounds');
ylabel('Alive nodes');
title('Nodes alive during rounds');

figure(4);
plot(r,P3.PACKETS_TO_BS3,'--k');
legend('I-SEP');
xlabel('Number of rounds');
ylabel('Throughput');
title('Packets sent to the base station');





