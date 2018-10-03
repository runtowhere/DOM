clear;


tmp1=[0.1838670,0.1838670,0.1838670,0.6950514,0.6950514,0.9656013];
ksi=[tmp1,-tmp1,-tmp1,tmp1];

tmp2=[0.1838670,0.6955014,0.9656013,0.1838670,0.6950514,0.1838670];
eta=[tmp2,tmp2,-tmp2,-tmp2];

tmp3=[0.1609517*2,0.3626469*2,0.1609517*2,0.3626469*2,0.3626469*2,0.1609517*2];
w=[tmp3,tmp3,tmp3,tmp3];

k=1.0;
ks=5.0;  %scattering coefficient
Lx=1;
Ly=1;
nx=40;
ny=40;
dx=Lx/nx;
dy=Ly/ny;
Tg=1000;
sigma=5.67e-8;

Ib=sigma*Tg^4/pi;

ndir=length(ksi);

n1=[1,0];
n2=[0,1];
n3=[-1,0];
n4=[0,-1];


dir=1;

tic;

Intensity(nx+2,ny+2,ndir)=0;
na=(nx+2)*(ny+2);

I=zeros(na,1);
A=zeros(na,na);
b=zeros(na,1);

G2=zeros(nx+2,ny+2);
G1=G2;

for it=1:30
    
    G3=G2;
    for id=1:ndir

        %firstly get the incoming intensity at the boundary;

    %     if ksi(id)>0 && eta(id)>0   %第一象限
    %         
    %     end
        
     
        
        for i=1:nx+2
            for j=1:ny+2
                li=i+(j-1)*(nx+2);
                if i==1 || i==nx+2 || j==1 || j==ny+2
                    A(li,li)=1.0;
                    b(li)=0;
                    continue;
                end

                %与象限无关
         %       A(li,li)=abs(ksi(id))*dy+abs(eta(id))*dx+k*dx*dy;

                li_w=li-1;  %左面和下面的辐射强度值的系数，在大矩阵A中都在同一行！
                li_s=li-(nx+2);
                li_e=li+1;
                li_n=li+(nx+2);

                aw=max(-ksi(id)*-1,0.0);
                ae=max(-ksi(id)*1.0,0.0);
                as=max(-eta(id)*-1,0.0);
                an=max(-eta(id)*1.0,0.0);

                A(li,li)=aw*dy+ae*dy+as*dx+an*dx+(k+ks)*dx*dy; 
                A(li,li_w)=-aw*dy;
                A(li,li_e)=-ae*dy;
                A(li,li_s)=-as*dx;
                A(li,li_n)=-an*dx;

                b(li)=k*Ib*dx*dy+ks/4/pi*G2(i,j)*dx*dy;

            end
        end

        I=A\b;   % solve the RTE in id-th direction 

        for i=2:nx+1
        for j=2:ny+1
            li=i+(j-1)*(nx+2);       
            Intensity(j,i,id)=I(li);
        end
        end
    end
    
    for i=2:nx+1
        for j=2:ny+1
            s=0;
            for m=1:ndir
                s=s+Intensity(i,j,m)*w(m);
            end
            G2(i,j)=s;   %update incident radiation
        end
    end
    relative_error=abs(G2-G3)./(G3+1.0e-10);
    max_error=max(max(relative_error));
    
    it,max_error
    if max_error<1.0e-7      
       break;
    end
    
    
end





x=linspace(dx/2,Lx-dx/2,nx);
y=linspace(dx/2,Lx-dx/2,nx);
[X,Y]=meshgrid(x,y);


for i=2:nx+1
    for j=1:ny+1
        for m=1:ndir
            G1(i,j)=G1(i,j)+Intensity(i,j,m)*w(m);
        end
    end
end


contourf(X,Y,G2(2:nx+1,2:nx+1),'ShowText','on');
axis equal;

toc;