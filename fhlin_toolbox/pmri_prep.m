function [A,profile_est,C,PRIOR]=pmri_prep(varargin)
%
%	pmri_prep		prepare for SENSE reconstruction
%
%
%	[A,S,C,P]=pmri_prep('ref',ref,'obs',obs,'sample_vector',sample_vector);
%
%	INPUT:
%	ref: input refernce data of [n_PE, n_FE, n_chan]. 
%		n_PE: # of phase encoding before acceleration
%		n_FE: # of frequency encoding
%		n_chan: # of channel
%	obs: input accelerated data of [n_PE, n_FE, n_chan]. 
%		n_PE: # of phase encoding after acceleration
%		n_FE: # of frequency encoding
%		n_chan: # of channel
%	sample_vector: 1D vector with entries of 0 or 1: [1, n_PE].
%		n_PE: # of phase encoding steps before acceleration
%		"0" indicates the correponding entries are not sampled in accelerated scan.
%		"1" indicates the correponding entries are sampled in accelerated scan.
%	'flag_display': value of either 0 or 1
%		It indicates of debugging information is on or off.
%	
%
%	OUTPUT:
%	A: 2D fourier aliasing matrix of [n_PE_acc, n_PE].
%		n_PE_acc: # of accelerated phase encoding steps
%		n_PE: # of phase encoding steps before acceleration
%	S: coil sensitivity maps of [n_PE, n_FE, n_chan].
%		n_PE: # of phase encoding
%		n_FE: # of frequency encoding
%		n_chan: # of channel
%	C: noise covariance matrix of [n_chan, n_chan].
%		n_chan: # of channel
%	P: prior of [n_PE, n_FE].
%		n_PE: # of phase encoding
%		n_FE: # of frequency encoding
%
%---------------------------------------------------------------------------------------
%	Fa-Hsuan Lin, Athinoula A. Martinos Center, Mass General Hospital
%
%	fhlin@nmr.mgh.harvard.edu
%
%	fhlin@jan. 20, 2005


C=[];
flag_ivs=0;

ref=[];
obs=[];
profile_est=[];

opt_level=[];

flag_profile_wavelet=0;
flag_profile_polynomial=1;

flag_display=0;
flag_sense_1d=1;

flag_profile_reim=0;
flag_profile_maph=1;


for i=1:floor(length(varargin)/2)
    option=varargin{i*2-1};
    option_value=varargin{i*2};
    switch lower(option)
	case 'ref'
		ref=option_value;
	case 'obs'
		obs=option_value;
	case 'sample_vector';
		sample_vector=option_value;
	case 'profile_opt'
		profile_opt=option_value;
	case 'opt_level'
		opt_level=option_value;
	case 'flag_ivs'
		flag_ivs=option_value;
	case 'flag_profile_wavelet'
		flag_profile_wavelet=option_value;
	case 'flag_profile_polynomial'
		flag_profile_polynomial=option_value;
	case 'flag_display'
		flag_display=option_value;
	case 'flag_sense_1d'
		flag_sense_1d=option_value;
	case 'flag_sense_2d'
		flag_sense_2d=option_value;
	case 'flag_profile_maph'
		flag_profile_maph=option_value;
	case 'flag_profile_reim'
		flag_profile_reim=option_value;
	otherwise
	        fprintf('unknown option [%s]!\n',option);
	        keyboard;
        	fprintf('error!\n');
        return;
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_ref=size(ref,3);
n_obs=size(obs,3);
if(n_ref<2 | n_obs <2)
	fprintf('Error! Reference/observation data is either empty or just single coil!\n');
	return;
end;
if(n_ref ~= n_obs)
	fprintf('Error! The number of channels in reference/observation data size is unequal!\n');
	return;
end;	

%the number of RF array coil
n_coil=n_obs; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	preparing refernce data: estimating coil sensitivity maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(flag_display)
	fprintf('preparing observation data...\n');
end;

rms_ref=zeros(size(ref(:,:,1)));
for i=1:n_coil
	rms_ref=rms_ref+abs(ref(:,:,i)).^2;
end;
rms_ref=sqrt(rms_ref./n_coil);
mask_ref = rms_ref > 0.1.*max(rms_ref(:));  %% a mask

if(flag_ivs)
	est=0;
else
	if(isempty(profile_est))
		est=1;
	else
		est=0;
	end;
end;

for i=1:n_coil
	if(est)
		%estimating B1 map using DWT/MVP
		if(flag_profile_wavelet)
       			%estimating profile by wavelet
			max_level=ceil(log10(max(size(ref,1)))/log10(2))-1;
			[ropt,popt,lopt]=correct111501(abs(ref(:,:,i)),'opt_level',opt_level);
			profile_est(:,:,i)=popt;
					
			profile_est(:,:,i)=fmri_scale(profile_est(:,:,i),1,0);
       
			%superimpose phase information
			profile_est(:,:,i)=profile_opt(:,:,i).*exp(sqrt(-1.0).*angle(ref(:,:,i)));

			%estimating B1 map using polynomial fitting
		elseif(flag_profile_polynomial)
			if(flag_profile_reim+flag_profile_maph==0)
				fprintf('error! must set either estimatinng magnitude/phase or real/imag profile!\n');
			end;
			if(flag_profile_reim+flag_profile_maph==2)
				fprintf('error! must set either estimatinng magnitude/phase or real/imag profile!\n');
			end;
			
			[dd,profile_est(:,:,i)]=pmri_poly_sensitivity(ref(:,:,i),3,'mask',mask_ref,'flag_reim',flag_profile_reim,'flag_maph',flag_profile_maph);
		end;
	else
		profile_est(:,:,i)=ref(:,:,i);
	end;

end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	preparing aliasing matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(flag_display)
	fprintf('preparing aliasing matrix...\n');
end;

if(flag_sense_1d)
	if(flag_display)
		fprintf('1D sense acceleration...\n');
	end;

	if(isempty(sample_vector))
		fprintf('Error! Must provide [sample_vector] variable to create aliasing matrix.\n');
		return;
	else
		A=pmri_alias_matrix(sample_vector);
	end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	preparing noise covariance matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(isempty(C))
	if(isempty(obs))	%no data
		fprintf('no data to build noise covariance matrix model!\nerror!\n');
		return;
	else
		%C=eye(size(obs,3));
		tmp=zeros(size(ref,1),size(ref,2));
		for i=1:size(ref,3)
			tmp=tmp+abs(ref(:,:,i).^2);
			tmp2=ref(:,:,i);
			yy(:,i)=tmp2(:);
		end;
		tmp=sqrt(tmp./size(ref,3));
		
		mask=zeros(size(tmp));
		noise_idx=find(tmp<mean(tmp(:))/10);
		mask(noise_idx)=1;
		
		noise_idx=[1:size(obs,1)*size(obs,2)];
		
		C=cov(yy(noise_idx,:));
	end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	preparing prior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(flag_display)
	fprintf('preparing prior...\n');
end;

%reference full-FOV image from all array elements: RMS image with no phase
PRIOR=zeros(size(ref(:,:,1)));
if(~flag_ivs)
	channel_power=squeeze(sum(sum(abs(ref).^2,1),2));
	for i=1:n_coil
		%PRIOR=PRIOR+ref(:,:,i)./profile_est(:,:,i);
		PRIOR=PRIOR+ref(:,:,i)./channel_power(i);
	end;
	PRIOR=PRIOR./n_coil;

	for i=1:n_coil
		xx(:,:,i)=A*(PRIOR.*profile_est(:,:,i));
	end;
	alpha=(xx(:)'*obs(:))/(xx(:)'*xx(:));
	PRIOR=PRIOR.*alpha;
else
	PRIOR=ones(size(ref(:,:,1)));	
end;

