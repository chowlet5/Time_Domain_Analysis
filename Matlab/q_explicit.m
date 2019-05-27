%**************************************************************************
%This function computes generalized coordinates using the Fourth Order
% Runge-Kutta explicit scheme
%**************************************************************************
function[q,qd,qdd] =q_explicit(M,P,C,K,k,N,dt)
for i=1:k
    for j=1:N
        if j==1
            q(i,j)=0.0;                              % set initial displacement value
			qd(i,j)=0.0;                             % set initial velocity value
    		qdd(i,j)=(1.0/M(i,i))*(P(i,j)-C(i,i)*qd(i,j)-K(i,i)*q(i,j));
        else
            q1=q(i,j-1);
			qd1=qd(i,j-1);
			qdd1=(1.0/M(i,i))*(P(i,j-1)-C(i,i)*qd1-K(i,i)*q1);
              
            q2=q(i,j-1) + 0.5*dt*qd1;
			qd2=qd(i,j-1)+ 0.5*dt*qdd1;
			qdd2=(1.0/M(i,i))*(P(i,j-1)-C(i,i)*qd2-K(i,i)*q2);
              
            q3=q(i,j-1) + 0.5*dt*qd2;
			qd3=qd(i,j-1) + 0.5*dt*qdd2;
			qdd3=(1.0/M(i,i))*(P(i,j-1)-C(i,i)*qd3-K(i,i)*q3);

            q4=q(i,j-1) + dt*qd3;
			qd4=qd(i,j-1) + dt*qdd3;
			qdd4=(1.0/M(i,i))*(P(i,j-1)-C(i,i)*qd4-K(i,i)*q4);
              
            q(i,j)=q(i,j-1) + (dt/6.0)*(qd1+ 2*qd2 + 2*qd3 + qd4);
			qd(i,j)=qd(i,j-1) + (dt/6.0)*(qdd1+ 2*qdd2 + 2*qdd3 + qdd4);
			qdd(i,j)=(1.0/M(i,i))*(P(i,j)-C(i,i)*qd(i,j)-K(i,i)*q(i,j));
            
        end
    end
end

