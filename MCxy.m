function MCxy()
%preparation
n=input('linear dimension');
beta=input('-J/T');
numoftrial=input('MCStep');
numofva=0;numofvb=0;
link={};
for i=0:(n-1)
    for j=0:(n-1)
        link{i+j*n+1}=[mod(j-1,n)*n+i,mod(i-1,n)+j*n,mod(j+1,n)*n+i,mod(i+1,n)+j*n]+1;
    end
end

%intialize
orien=zeros(1,n^2);
orien=unifrnd(0,2*pi,1,n^2);

%main program

for num=1:numoftrial
    numofva=0;numofvb=0;
    loopupdate;
    
    set(rectangle('position',[0.5,0.5,n+1,n+1]),'EdgeColor','black')
    title('XY Model')
    xlabel('x')
    ylabel('y')
    hold on
    for i=0:(n-1)
        for j=0:(n-1)
            plot(i+1,j+1,'.')
            %hold on
            draw_arrow([i+1,j+1],0.5*[cos(orien(i+j*n+1)),sin(orien(i+j*n+1))]+[i+1,j+1],2.0)
            %hold on
            v=identifyV(i,j);
            draw_box(i+1,j+1,v)            
            %hold on
        end
    end
    text(n+1-n/10,n+1+n/10,['NumofVortice:',num2str(numofva);'NumofAntiVor:',num2str(numofvb)])
    text(n/2,n+n/10+n/50,'by Yuliang Huang')
    hold off
    drawnow
    clf
    num2str([num,numofva,numofvb])
    
end
%MCstep
    function loopupdate()
        cluster=[];
        stack=[];
        theta=unifrnd(0,2*pi,1);
        x=unidrnd(n,1)-1;y=unidrnd(n,1)-1;
        cluster(1)=x+y*n+1;
        stack(1)=x+n*y+1;
        nstack=1;
        ncluster=1;
        orien(x+y*n+1)=mod(2*theta+pi-orien(x+n*y+1),2*pi);
        while nstack~=0
                ia=unidrnd(nstack,1);
                t=stack(ia);
                for ib=1:4
                    t1=link{t}(ib);
                    if ismember(t1,cluster)==0
                        pr=1-exp(min(0,2*beta*cos(orien(t)-theta)*cos(orien(t1)-theta)));
                        if unifrnd(0,1,1)<pr
                            orien(t1)=mod(2*theta+pi-orien(t1),2*pi);
                            cluster(ncluster+1)=t1;
                            ncluster=ncluster+1;
                            stack(nstack+1)=t1;
                            nstack=nstack+1;
                        end
                    end
                end
                stack(ia)=[];
                nstack=nstack-1;
        end
    end
    
%submodule for plot
    function draw_arrow(startpoint,endpoint,headsize)
        %by Ryan Molecke 
        % accepts two [x y] coords and one double headsize 
        v1 = headsize*(startpoint-endpoint)/2.5; 
        theta = 22.5*pi/180; 
        theta1 = -1*22.5*pi/180; 
        rotMatrix = [cos(theta)  -sin(theta) ; sin(theta)  cos(theta)];
        rotMatrix1 = [cos(theta1)  -sin(theta1) ; sin(theta1)  cos(theta1)];  
        v2 = v1*rotMatrix; 
        v3 = v1*rotMatrix1; 
        x1 = endpoint;
        x2 = x1 + 0.5*v2; 
        x3 = x1 + 0.5*v3; 
        hold on; 
        fill([x1(1) x2(1) x3(1)],[x1(2) x2(2) x3(2)],[0 0 0]);% this fills the arrowhead (black) 
        plot([startpoint(1) endpoint(1)],[startpoint(2) endpoint(2)],'linewidth',2,'color',[0 0 0]);
    end
    function b=identifyV(x,y)
        delta=0;
        angles=orien([x+n*y,mod(x+1,n)+y*n,mod(x+1,n)+mod(y+1,n)*n,x+mod(y+1,n)*n]+1);
        for ic=1:3
            temp=angles(ic+1)-angles(ic);
            if abs(temp)<pi
                delta=delta+temp;
            elseif temp>pi
                delta=delta+temp-2*pi;
            else
                delta=delta+temp+2*pi;
            end
        end
        if delta>pi
            b=1;
        elseif delta<-pi
            b=-1;
        else
            b=0;
        end
    end
    function draw_box(x,y,v)
        if v==1
            plot([x,x+0.5],[y,y],'r')
            hold on
            plot([x,x],[y,y+0.5],'r')
            plot([mod(x,n)+1,mod(x,n)+0.5],[y,y],'r')
            plot([mod(x,n)+1,mod(x,n)+1],[y,y+0.5],'r')
            plot([x,x],[mod(y,n)+1,mod(y,n)+0.5],'r')
            plot([x,x+0.5],[mod(y,n)+1,mod(y,n)+1],'r')
            plot([mod(x,n)+1,mod(x,n)+0.5],[mod(y,n)+1,mod(y,n)+1],'r')
            plot([mod(x,n)+1,mod(x,n)+1],[mod(y,n)+1,mod(y,n)+0.5],'r')
            numofva=numofva+1;
        elseif v==-1
            plot([x,x+0.5],[y,y],'b')
            hold on
            plot([x,x],[y,y+0.5],'b')
            plot([mod(x,n)+1,mod(x,n)+0.5],[y,y],'b')
            plot([mod(x,n)+1,mod(x,n)+1],[y,y+0.5],'b')
            plot([x,x],[mod(y,n)+1,mod(y,n)+0.5],'b')
            plot([x,x+0.5],[mod(y,n)+1,mod(y,n)+1],'b')
            plot([mod(x,n)+1,mod(x,n)+0.5],[mod(y,n)+1,mod(y,n)+1],'b')
            plot([mod(x,n)+1,mod(x,n)+1],[mod(y,n)+1,mod(y,n)+0.5],'b')
            numofvb=numofvb+1;
        else
            return
        end
    end

end