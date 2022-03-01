function res=My_comple(S,G,truth,c,omega,b,p,mode,r,lambda)
%% Initialization
[~,n2,~]=size(S);
for i=1:n2
    S_sss(:,:,i)=G{i}'*S{i}*G{i};
    Q2{i}=zeros(size(S{i}));
    P{i}=Q2{i};
    E{i}=P{i};
end
dim=size(S_sss);
[n1,n2,n3]=size(S_sss);
mu=1e-4;
tol = 1e-8; max_iter = 500; rho = 1.1;

Q1=zeros(dim); Q3=zeros(dim); W=zeros(dim); Y=zeros(dim);
Z=zeros(dim); M=zeros(dim); F=zeros(n2,c); U=zeros(n2,c);
iter=0;

alpha=1/n3*ones(1,n3);

Z=S_sss;
beta=ones(n1,1);

for i=1:n3
    [nu,~]=size(S{i});
    omega(:,:,i)=G{i}'*ones(nu,nu)*G{i};%ones(numFold,numFold);
end
for i=1:max_iter
    iter=iter+1;
    X_k=Z;
    Z_k=M;
    for j=1:n3
        P{j}=S{j}+1/mu*Q2{j}-E{j};
    end
    %% Update Z
    for j=1:n3
        Z(:,:,j)=0.5*(M(:,:,j)-1/mu*Q1(:,:,j)+W(:,:,j)-1/mu*Q3(:,:,j))-1/mu*(Y(:,:,j).*omega(:,:,j));
    end
    
    %% Update  W
    for j=1:n3
        temp1 = L2_distance_1(F',F');
        temp2 = Z(:,:,j)+Q3(:,:,j)/mu;
        linshi_W = temp2-alpha(j)^r*b*temp1/mu;
        linshi_W = linshi_W-diag(diag(linshi_W));
        for ic = 1:size(Z(:,:,j),2)
            ind = 1:size(Z(:,:,j),2);
            %             ind(ic) = [];
            W(ic,ind,j) = EProjSimplex_new(linshi_W(ic,ind));
        end
    end
    clear temp1 temp2
    
    %% Update M
    [M,~,~] = prox_tnn(Z+Q1/mu,beta/mu,p,mode);
    
    %% Update F
    temp_W=zeros(size(W(:,:,1)));
    for j=1:n3
        temp_W=alpha(j)^r*W(:,:,j)+temp_W;
    end
    temp_W = (temp_W+temp_W')/2;                                                        %update F
    L_D = diag(sum(temp_W));
    L_Z = L_D-temp_W;
    F_old = F;
    [F, ~, ~]=eig1(L_Z, c, 0);
    
    %% Update E
    for j=1:n3
        temp1 = S{j}-P{j}+Q2{j}/mu;
        temp2 = lambda/mu;
        E{j}= max(0,temp1-temp2)+min(0,temp1+temp2);
    end
    
    clear temp1 temp2
    %% Update alpha
    for j=1:n3
        h(j)=trace(F'*W(:,:,j)*F);
        temp(j)=((r*h(j))^(1/(1-r)));
    end
    if sum(temp)==0
        for j=1:n3
            alpha(j)=1/n3;
        end
    elseif sum(temp)~=0
        for j=1:n3
            alpha(j)=temp(j)/sum(temp);
        end
    end
    %% Checking Convergence
    chgX=max(abs(Z(:)-X_k(:)));
    chgZ=max(abs(M(:)-Z_k(:)));
    chgX_Z=max(abs(Z(:)-M(:)));
    chg=max([chgX chgZ chgX_Z]);
    
    if iter == 1 || mod(iter, 10) == 0
        disp(['iter ' num2str(iter) ', mu = ' num2str(mu) ', chg = ' num2str(chg) ', chgX = ' num2str(chgX) ', chgZ = ' num2str(chgZ) ',chgX_Z = ' num2str(chgX_Z) ]);
    end
    
    if chg<tol
        break;
    end
    %% Update Lagrange multiplier
    for j=1:n3
        tP(:,:,j)=G{j}'*P{j}*G{j};
    end
    Q1=Q1+mu*(Z-M);
    for j=1:n3
        Q2{j}=Q2{j}+mu*(S{j}-P{j}-E{j});
    end
    Q3=Q3+mu*(Z-W);
    Y=Y+mu*(Z-tP).*omega;
    mu=rho*mu;

    %% Clustering
    new_F = F;
    norm_mat = repmat(sqrt(sum(new_F.*new_F,2)),1,size(new_F,2));
    for i = 1:size(norm_mat,1)
        if (norm_mat(i,1)==0)
            norm_mat(i,:) = 1;
        end
    end
    new_F = new_F./norm_mat;
    repeat = 5;
    for iter_c = 1:repeat
        pre_labels    = kmeans(real(new_F),c,'emptyaction','singleton','replicates',20,'display','off');
        result_LatLRR = ClusteringMeasure(truth, pre_labels);
        AC(iter_c)    = result_LatLRR(1)*100;
        MIhat(iter_c) = result_LatLRR(2)*100;
        Purity(iter_c)= result_LatLRR(3)*100;
    end
    mean_ACC = mean(AC);
    mean_NMI = mean(MIhat);
    mean_PUR = mean(Purity);
    res=[mean_ACC mean_NMI mean_PUR];
    if iter == 1 || mod(iter, 10) == 0
        disp(['iter ' num2str(iter) ', mean_ACC = ' num2str(mean_ACC) ', mean_NMI = ' num2str(mean_NMI) ', mean_PUR = ' num2str(mean_PUR)]);
        disp(['---------------------------------------------------------------------------------------'])
    end
end
