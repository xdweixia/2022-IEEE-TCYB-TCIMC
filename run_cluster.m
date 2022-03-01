clear;
clc;

addpath(genpath('funs/'));
addpath(genpath('tSVD/'));

res = [];
[res_IMG]=[];
[res_PVC]=[];
[res_GPVCO]=[];

%% Load ORL dataset
f=1;
load('dataset\ORL_4views.mat'); c=40;truth=truth';
load('dataset\ORL_e2_fold.mat');
ind_folds = folds{f};
numClust = length(unique(truth));
num_view = length(X);
[numFold,numInst]=size(ind_folds);
fid=fopen('ORL_Results.txt','a');


result=[];

for iv = 1:num_view
    X1 = X{iv}';
    X1 = NormalizeFea(X1,1);
    ind_0 = find(ind_folds(:,iv) == 0);
    X1(ind_0,:) = [];
    Y{iv} = X1';
    W1 = eye(size(ind_folds,1));
    W1(ind_0,:) = [];
    G{iv} = W1;
end
%% Graph construction
S_temp=graph_construction(Y);
for i=1:num_view
    S(:,:,i)=G{i}'*S_temp{i}*G{i};
    [nu,~]=size(S_temp{i});
    omega(:,:,i)=G{i}'*ones(nu,nu)*G{i};
end

%% Training
for r = 1.3
    for p= 0.6
        for beta = 0.3
            mode=2;
            [res] = My_comple(S_temp,G,truth,c,omega,beta,p,mode,r,0);
            [result]=[result;res];
            fprintf(fid,'r: %f ',r);
            fprintf(fid,'p: %f ',p);
            fprintf(fid,'beta: %f ',beta);
            fprintf(fid,'res1: %g %g %g  \n',res(:,1:3));
        end
    end
end