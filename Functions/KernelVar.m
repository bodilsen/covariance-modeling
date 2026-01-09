function [V,bias] = KernelVar(H,n,m,IQ,IV,om2,Vers,topF)
    
	w = MKerWeight(H, Vers, topF)';
    w(1) = w(1)+ 0.5;

	V    = KernelVariance(w,n,IQ,IV,om2,2,m);
    bias = (w(1)-w(2))*2*om2*n;
end


function [V] = KernelVariance(w,n,IQ,IV,ome_sq,lamsq,m)

	V = 4*(IQ*KA(w)/n+2*ome_sq*IV*KB(w)+ome_sq^2*(n*KC(w,lamsq)+KD(w,lamsq)+KE(w,lamsq,m)/m));
end


function [A] = MatrixA(H)

    A = toeplitz([1 zeros(1,H)]); A(1,1) = 1/2;
end

function [B] = MatrixB(H)

    B = toeplitz([2 -1 zeros(1,H-1)]); B(1,1) = 1;
end

function [C] = MatrixC(H)
% The C-matrix (when lambda^2=0)
    C  = toeplitz([6 -4 1 zeros(1,H-2)]);
	C(1:2,1:2)=[1 -2; -2 5];
end

function [D] = MatrixD(H)
% The D-matrix (when lambda^2=0)
	H2=H+1; 	
    H=max(H,3);
    D=zeros(H+1,H+1);
	D(1:3,1:3)= [ -2, 4  ,-2; 4, -21/2, 9; -2, 9, -15];
	D(4,2)=-5/2;
	D(2,4)=-5/2;
	D(4,3)=11;
	D(3,4)=11;

	for i=4:H+1
		D(i,i)   = D(i-1,i-1)-3;
		D(i-1,i) = D(i-2,i-1)+2;
		D(i,i-1) = D(i-2,i-1)+2;
		D(i-2,i) = D(i-3,i-1)-1/2;
		D(i,i-2) = D(i-3,i-1)-1/2;
	end
	D =  D(1:H2,1:H2);
end

function [E] = MatrixE(H)
%  The E-matrix (when lambda^2=0 and m=1)
    E  = toeplitz([8 -5 1 zeros(1,H-2)]);
	E(1:2,1:2) = [2 -3; -3 7];
end
   

function [KA] = KA(w)
% computes w'Aw
    if size(w,2)>1, w=w';end; 
	H = length(w)-1;
	w = H^(-1/2).*w;	% terms divided by H for numerical reasons
	S = w'*w-w(1)^2/2;
	KA = S*H;
end

function [KB] = KB(w)
% computes w'Bw
    if size(w,2)>1, w=w'; end; 
	H = length(w)-1;
	w = H^(1/2).*w;     % terms multiplied by H for numerical reasons
	if (H<5)
		S = w'*MatrixB(H)*w;
    else	
		S = w(1)*(w(1)-w(2));
		for h=2:H
			S = S-w(h)*(w(h-1)-2*w(h)+w(h+1));
        end
		S = S+w(H+1)*(-w(H)+2*w(H+1));
    end
	KB = S/H;
end

function [KC] = KC(w,lamsq)	  	 
% computes w'Cw
    if size(w,2)>1, w=w'; end; 
	H = length(w)-1;
	w=H.*w;	 
	if (H<5)
		S = w'*MatrixC(H)*w;
    else	
		S = w(1)*(1*w(1)-2*w(2)+w(3));
		S = S+w(2)*(-2*w(1)+5*w(2)-4*w(3)+w(4));
		for h=3:H-1
			S = S+w(h)*(w(h-2)-4*w(h-1)+6*w(h)-4*w(h+1)+w(h+2));
        end
		S = S+w(H)*(w(H-2)-4*w(H-1)+6*w(H)-4*w(H+1));
		S = S+w(H+1)*(w(H-1)-4*w(H)+6*w(H+1));
    end
	S  = S+lamsq*(w(1)-w(2))^2;
    KC = S/H/H;
end

function [KD] = KD(w,lamsq)
% computes w'Dw
    if size(w,2)>1, w=w'; end; 
	H = length(w)-1;
	w = sqrt(H).*w;     % terms multiplied by H for numerical reasons
	if (H<5)
		S=w'*MatrixD(H)*w;
    else
		S=w(1)*(-2*w(1)+4*w(2)-2*w(3));
		S=S+w(2)*(4*w(1)-10.5*w(2)+9*w(3)-5/2*w(4));
		for h=3:H-1
			S=S+w(h)*(-(h+1)/2*w(h-2)+(2*(h-1)+5)*w(h-1)-(3*(h-1)+9)*w(h)+(2*(h-1)+7)*w(h+1)-(h+3)/2*w(h+2));
        end
		S = S+w(H)*(-(H-1+2)/2*w(H-2)+(2*(H-1)+5)*w(H-1)-(3*(H-1)+9)*w(H)+(2*(H-1)+7)*w(H+1));
		S = S+w(H+1)*(-(H+2)/2*w(H-1)+(2*H+5)*w(H)-3*(H+3)*w(H+1));
    end
	S  = S-lamsq*(w(1)-w(2))^2;
	KD = S/H;
end

%fprintf('%10.6g\n',S);

function [KE] = KE(w,lamsq,m)
% computes w'Cw
    if size(w,2)>1, w=w'; end; 
	H = length(w)-1;
	w = sqrt(H).*w;     % terms multiplied by H for numerical reasons
    
	if (H<5)
		S = w'*MatrixE(H)*w;
    else
		S = w(1)*(2*w(1)-3*w(2)+w(3));
		S = S+w(2)*(-3*w(1)+7*w(2)-5*w(3)+w(4));
		for h=3:H-1
			S=S+w(h)*(w(h-2)-5*w(h-1)+8*w(h)-5*w(h+1)+w(h+2));
        end
		S = S+w(H)*(w(H-2)-5*w(H-1)+8*w(H)-5*w(H+1));
		S = S+w(H+1)*(w(H-1)-5*w(H)+8*w(H+1));
    end
	S  = S+(lamsq/2+m-1)/m/m*(w(1)-w(2))^2;
	KE = S/H;
end
    

