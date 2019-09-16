function [x_unreg, g_unreg, x_reg, g_reg]=pmri_core(varargin);
%
%	pmri_core		perform SENSE reconstruction
%
%
%	[recon_unreg, g_unreg]=pmri_core('A',A,'Y',Y,'S',S,'C',C,'flag_unreg',1,'flag_unreg_g',1);
%
%	INPUT:
%	A: 2D fourier aliasing matrix of [n_PE_acc, n_PE].
%		n_PE_acc: # of accelerated phase encoding steps
%		n_PE: # of phase encoding steps before acceleration
%	Y: input data of [n_PE, n_FE, n_chan]. 
%		n_PE: # of phase encoding
%		n_FE: # of frequency encoding
%		n_chan: # of channel
%	S: coil sensitivity maps of [n_PE, n_FE, n_chan].
%		n_PE: # of phase encoding
%		n_FE: # of frequency encoding
%		n_chan: # of channel
%	C: noise covariance matrix of [n_chan, n_chan].
%		n_chan: # of channel
%	P: prior of [n_PE, n_FE].
%		n_PE: # of phase encoding
%		n_FE: # of frequency encoding
%	'flag_unreg': value of either 0 or 1
%		reconstruct un-regularized SENSE image
%	'flag_reg': value of either 0 or 1
%		reconstruct regularized SENSE image
%	'flag_unreg_g': value of either 0 or 1
%		computer un-regularized g-factor map
%	'flag_reg_g': value of either 0 or 1
%		computer regularized g-factor map
%	'flag_display': value of either 0 or 1
%		It indicates of debugging information is on or off.
%
%
%	OUTPUT:
%	recon_unreg: 2D un-regularized SENSE reconstruction [n_PE, n_PE].
%		n_PE: # of phase encoding steps
%		n_FE: # of frequency encoding steps
%	g_unreg: 2D un-regularized g-factor map [n_PE, n_PE].
%		n_PE: # of phase encoding steps
%		n_FE: # of frequency encoding steps
%
%---------------------------------------------------------------------------------------
%	Fa-Hsuan Lin, Athinoula A. Martinos Center, Mass General Hospital
%
%	fhlin@nmr.mgh.harvard.edu
%
%	fhlin@jan. 20, 2005

A=[];

S=[];

C=[];

Y=[];

P=[];

sample_vector=[];

flag_display=0;

flag_reg=0;
flag_reg_g=0;

flag_unreg=1;
flag_unreg_g=0;

flag_lsqr=0;

for i=1:floor(length(varargin)/2)
    option=varargin{i*2-1};
    option_value=varargin{i*2};
    switch lower(option)
    case 'a'
        A=option_value;
    case 's'
        S=option_value;
    case 'c'
        C=option_value;
    case 'p'
        P=option_value;
    case 'y'
        Y=option_value;
    case 'lambda'
    	lambda=option_value;
    case 'flag_display'
	flag_display=option_value;
    case 'reg_lambda'
        reg_lambda=option_value;
    case 'flag_reg'
	flag_reg=option_value;
    case 'flag_reg_g'
	flag_reg_g=option_value;
    case 'flag_unreg'
	flag_unreg=option_value;
    case 'flag_lsqr'
	flag_lsqr=option_value;
    case 'flag_unreg_g'
	flag_unreg_g=option_value;
    case 'recon_channel_weight'
	recon_channel_weight=option_value;
    otherwise
        fprintf('unknown option [%s]!\n',option);
        fprintf('error!\n');
        return;
    end;
end;


x_unreg=zeros(size(A,2),size(Y,2));
x_reg=zeros(size(A,2),size(Y,2));

g_unreg=zeros(size(A,2),size(Y,2));
g_reg=zeros(size(A,2),size(Y,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	preparation of noise covariance and whitening matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(isempty(C))
	if(isempty(Y))	%no data
		fprintf('no data to build noise covariance matrix model!\nerror!\n');
		return;
	else
		C=eye(size(Y,3));
	end;
end;
%prepare whitening data
[u,s,v]=svd(C);
W=pinv(sqrt(s))*u';	%whitening matrix

%ss=repmat(transpose(diag(C)),[size(Y,1),1]);
%N=diag(ss(:));
%[u_n,s_n,v_n]=svd(N);		
%W_n=sqrt(pinv(s_n))*u_n';


tic;

n_coil=size(Y,3);
Y=permute(Y,[1,3,2]);
S=permute(S,[1,3,2]);

for fe_idx=1:size(Y,3)
%for fe_idx=64:64
	if(flag_display)
		fprintf('.');
	end;
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	preparation of observation data
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	Y_w=Y(:,:,fe_idx)*W';
	
	Z=[];	
	for i=1:n_coil
		Z=cat(1,Z,Y_w(:,i));
	end;
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	preparation of encoding matrix
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	S_w=S(:,:,fe_idx)*W';
	
	E=[];
	for i=1:n_coil
		E=cat(1,E,A.*(ones(size(A,1),1)*transpose(S_w(:,i))));		
	end;
		
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	SENSE reconstruction using least squares
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if(flag_unreg)
		if(flag_lsqr)
			x_unreg(:,fe_idx)=lsqr(E,Z,0.1,size(E,2));
		else
			M=E'*E;
			x_unreg(:,fe_idx)=inv(M)*E'*Z;
		end;
	end;
	
	if(flag_unreg_g)
		M=E'*E;
		g_unreg(:,fe_idx)=diag(inv(M)).*diag(M);
	end;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	SENSE reconstruction using regularized least squares
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if(flag_reg|flag_reg_g)
		[u_e,s_e,v_e]=svd(E,0);
		
		Gamma=diag(diag(s_e)./(diag(s_e).^2+lambda(fe_idx).^2));
		
		Phi=diag(lambda(fe_idx).^2./(diag(s_e).^2+lambda(fe_idx).^2));
		
		if(flag_reg)
			x_reg(:,fe_idx)=v_e*Gamma*u_e'*Z+v_e*Phi*v_e'*P(:,fe_idx);
		end;
		
		if(flag_reg_g)
			g_reg(:,fe_idx)=diag(v_e*(Gamma).^2*v_e').*diag(v_e*(s_e).^2*v_e');
		else
			g_reg(:,fe_idx)=0;
		end;
	else
		x_reg(:,fe_idx)=0;
		g_reg(:,fe_idx)=0;
	end;
	
end;

if(flag_display)
	fprintf('\n');
	tt=toc;
	fprintf('SENSE recon. time=%2.2f\n',tt);
end;

