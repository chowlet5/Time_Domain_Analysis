%**************************************************************************
%This function computes generalized coordinates using the Newmark-beta    *
%linear acceleration implicit scheme                                      *
%**************************************************************************
function[q,qd,qdd] =q_implicit(M,P,C,K,k,N,dt)
ga = 0.5;                              %Gama for linear acceleration method
be = 0.1666667;                        %Beta for linear acceleration method

dpp=0.0;
dq=0.0;
dqd=0.0;
dqdd=0.0;

q = zeros(k,N);
qd= zeros(k,N);
qdd = zeros(k,N);
% Evaluation of constant intermediate matrix KP
for i=1:k
    for j=1:k
        KP(i,j)=K(i,j)+(ga/(be*dt))*C(i,j)+(1/(be*(dt*dt)))*M(i,j);
    end
end
% Evaluation of constant intermeiate matrix a
for i=1:k
    for j=1:k
        a(i,j)=(1/(be*dt))*M(i,j)+(ga/be)*C(i,j);
    end
end  
%Evaluation of constant intermediate matrix b
for i=1:k
    for j=1:k
        b(i,j)=(1/(2*be))*M(i,j)+dt*((ga/(2*be))-1)*C(i,j);
    end
end
% Repetitive calculations for step i
for i=1:k
    for j=1:N
        if j==1
            dpp=P(i,j)+a(i,i)*qd(i,j)+b(i,i)*qdd(i,j);
            
            dq=(1/KP(i,i))*dpp;
			dqd=(ga/(be*dt))*dq-(ga/be)*qd(i,j)+dt*(1-(ga/(2*be)))*qdd(i,j);
			dqdd=(1/(be*dt*dt))*dq-(1/(be*dt))*qd(i,j)-(1/(2*be))*qdd(i,j);
            
            q(i,j)=q(i,j)+dq;
            qd(i,j)=qd(i,j)+dqd;
            qdd(i,j)=qdd(i,j)+dqdd;
        else
            q(i,j)=q(i,j-1);
			qd(i,j)=qd(i,j-1);
			qdd(i,j)=qdd(i,j-1);
			
            dpp=P(i,j)-P(i,j-1)+a(i,i)*qd(i,j)+b(i,i)*qdd(i,j);
			dq=(1/KP(i,i))*dpp;
			dqd=(ga/(be*dt))*dq-(ga/be)*qd(i,j)+dt*(1-(ga/(2*be)))*qdd(i,j);
			dqdd=(1/(be*dt*dt))*dq-(1/(be*dt))*qd(i,j)-(1/(2*be))*qdd(i,j);
		
			q(i,j)=q(i,j-1)+dq;
			qd(i,j)=qd(i,j-1)+dqd;
			qdd(i,j)=qdd(i,j-1)+dqdd;
        end
    end
end
