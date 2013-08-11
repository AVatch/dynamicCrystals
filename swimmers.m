
clear all

%1s = tp*steps
steps=200;
deltat=1/steps;
runtime=10;

threshold=3;
N=5; %%number swimmers
V0_int=1;
 
for tp=[5]
    fileTitle=strcat('rho10__N600_tp',int2str(tp),'_rt250_V1_t1.avi');
    qq = 10;
    rho=qq/100;     %percent surface coverage
    
    %N=floor((4*box^2*rho)/(pi))
    %rho=pi*N/(4*box^2)
    box=floor(sqrt(N*pi/(4*rho)))
    %%%INITIAL
    part=box*rand(N,2);
    duration=runtime*tp/deltat;
    
    fprintf('Number of Particles:\t%d\n',N);
    fprintf('Box Size:\t\t\t\t%d\n',box);
    fprintf('Rho:\t\t\t\t\t%f\n',rho);
    fprintf('\t================\n')
    fprintf('\t================\n')
    
    diff_theta=0.01;%% initial conditions
    x=0*ones(N,duration);
    y=0*ones(N,duration);
    clear X;
    clear Y;
    coeff=0.87/2;
    thresh_contact=5;
    
    tic
    bb=0;

    for i=1:N
        theta(i)=2*pi*(rand()+0.5);
        x(i,1)=box*(rand());
        y(i,1)=box*(rand());
        V0(i)=V0_int;
    end

    k=1;

    for r=2:duration
        k=k+1;
        i=0;
        for l=1:N
            i=i+1;     
            x(i,k)=x(i,k-1)+deltat*cos(theta(i))*V0(i);
            y(i,k)=y(i,k-1)+deltat*sin(theta(i))*V0(i);
            V0(i)=V0_int;
            j=0;
            contact=0;
            for p=1:N
                j=j+1;
                dist(i,j)=sqrt((x(i,k)-x(j,k))^2+(y(i,k)-y(j,k))^2);
                if dist(i,j)>0
                    if dist(i,j)<threshold
                        tan_ang=(y(j,k)-y(i,k))/(x(j,k)-x(i,k));
         
                        if (x(j,k)-x(i,k))>0
                            ang=atan(tan_ang);
                        else
                            ang=atan(tan_ang)+pi;
                        end
            x(j,k)=x(j,k)-deltat*coeff/((dist(i,j)))^2*cos(ang);
            y(j,k)=y(j,k)-deltat*coeff/((dist(i,j)))^2*sin(ang);
            x(i,k)=x(i,k)+deltat*coeff/((dist(i,j)))^2*cos(ang);
            y(i,k)=y(i,k)+deltat*coeff/((dist(i,j)))^2*sin(ang);
         
            dist(i,j)=sqrt((x(i,k)-x(j,k))^2+(y(i,k)-y(j,k))^2);
%
%  	CORRECT FOR OVERLAP
%
            if dist(i,j)<1
                x(i,k)=x(i,k)-(1-dist(i,j))*cos(ang)/2;
                y(i,k)=y(i,k)-(1-dist(i,j))*sin(ang)/2;
                x(j,k)=x(j,k)+(1-dist(i,j))*cos(ang)/2;
                y(j,k)=y(j,k)+(1-dist(i,j))*sin(ang)/2;
            end
          end
       end
    end
theta(i)=theta(i)+sqrt(2/((steps)*tp))*(randn);

end


%%%%%%%PBC%%%%%
for j=1:N
    if x(j,k)>=box
        x(j,k)=x(j,k)-box;
    end
    if y(j,k)>=box
        y(j,k)=y(j,k)-box;
    end
    if y(j,k)<=0
        y(j,k)=y(j,k)+box;
    end
    if x(j,k)<=0
        x(j,k)=x(j,k)+box;
    end
end

if(mod(r,100)==0)
save(strcat(['test_',num2str(r),'frame'],'.txt'), 'x','y','-ascii')
end
    end
end

fprintf('\t================\n')

% aviobj = avifile(fileTitle,'compression','None');
% l=0;
%  for k=2:steps*tp/5:duration
%  l=l+1;
%  
%  
%  for i=1:N
%     %Plots particles
%      DrawCircle(x(i,k), y(i,k), 0.5,64, 'k');
%  
%     hold on;
%  end
%  
% 
%  hold off;
%  q=k/steps/tp;
%  title(q);
%  axis([0 box 0 box])
%  pause(0.1)
%  M=getframe;
%  %Get the image and put it into the movie
%  aviobj = addframe(aviobj,M)
%  end
%  aviobj = close(aviobj);
%  fprintf('\tSaved\n')
%  fprintf('\t================\n')
% 
% end
