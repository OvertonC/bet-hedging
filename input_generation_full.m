function input_generation_full
D=[];T=[];
Nsize = 8;
Njobs = 100000;
for i=1:1000000
    name = ['matrix',num2str(i)];
M = make_random_graph_strict(Nsize,4);
deg = sum(M,2);
D.(name) = var(deg);
T.(name) = M;
list(i) = var(deg);
mat(:,:,i) = M;
end

index = unique(list);
index2 = sort(index);
diff = find(index==index2);
for i = 1:length(index)
    j = find(list==index(i),1,'first');
    name = ['matrix',num2str(i)];
    out.(name) = mat(:,:,j);
end

for i = 1:length(index)
    name = ['matrix',num2str(i)];
    M = out.(name);
    degvar(i) = var(sum(M,2));
end

for i=1:Njobs
    out.seed=i;
    name = ['input',num2str(i-1),'.mat'];
    save(name,'out');
end

function [ G ] = make_random_graph_strict( N,n )


G=rand(N)<n/N;
G=tril(G);
G=G+G';
G=G-diag(diag(G));
while   sum(sum(G))~=N*n | find(sum(G)==0)~=0 | find(sum(G^2,2)==1)~=0 
    G=rand(N)<n/N;
    G=tril(G);
    G=G+G';
    G=G-diag(diag(G));
end
end

end

