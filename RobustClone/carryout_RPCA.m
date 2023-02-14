%clear
function [] = carryout_RPCA(input_file, output_file)
	addpath(genpath('matlab_and_R_scripts')); 
	tic
	%D=csvread('exampledata.csv',1,1);
	D=csvread(input_file,1,1);
	%% there are no unobservable entries
	%[A1,E1]= RPCA(D); % RPCA model

	%% there exist unobservable entries
	[m,n]=size(D); % In SNV data, 3 represents missing in general.
%                           In CNV data, the missing copy number at the particular locations can be replaced by a number that does not appear in the CNV matrix.  For example, the missing locations can be filled with -1 in D.
	ms=3; % ms represents missing data. In SNV data, if 3 represents missing, then ms=3; In CNV data, if -1 represents missing, then ms=-1.
	omega=find(D~=ms); 
	omegaC=find(D==ms);
	lambda=1/sqrt(max(m,n))*(1+3*length(omegaC)/(m*n));
	[A1,E1]= extendedRPCA(D,omega,lambda); % the extended RPCA model

	%% Integralization
	AA1=int8(A1);
	%save('example_RobustClone.mat','AA1')
	save(output_file,'AA1')
	toc
end
%function data=readdat(filename)
%	arguments
%      		filename {mustBeText}
%   	end
%end
%exit
