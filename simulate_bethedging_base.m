%%%% Main function_________________________________________________________
function simulate_bethedging_base
%%% This code provides the basics for simulating evolutionary graph theory
%%% in the presence of variable fitness. Users can specify contact
%%% networks, mean fitness, fitness variation, and initial conditions.
%%% Additionally, within the SimulationDBBV function, there is an option
%%% "var_type". If this is set to "Demographic", the model will use within
%%% generational variation. If set to "Environmental", the model will use
%%% between generational variation. At the end we provide examples of
%%% random graph generation functions, so that users can test different
%%% types of network. 


%%% Define mean fitness values
a=0.95; % Mean bet-hedger fitness
c=1; % Mean normal-type fitness


%%% Specify number of simulations to be run. Better to use low number and
%%% run in parallel.
tot = 10;

    
%%% Generate contact network
T = make_constant_n_symmetric_matrix(50,4);

%%% Specify scale parameters for gamma distributed fitness. These should be
%%% set up such that variance = mean*scale, (e.g. varb = variance_b/mean_b)
vara = 2.1053;
varb = 2.5;       

N = size(T,2); 
S=zeros(N,N);
R=zeros(N,N);

for I=1:tot
    % Specify initial conditions, "Ap" - bet-hedger location vector, "Bp" - 
    % normal-type location vector. Currently set up to generate random
    % vector with a single bet-hedger, but can change to any vector of
    % interest.
    Ap=make_random_vector(N,1); Bp=1-Ap; 
    
    [R]=SimulationDBBV(T,Ap,Bp,a,c,varb,vara); % Simulate the dynamics until fixation
    
    %%%% Aggregate output across multiple simulations______________________
    SR=length(R(1,:)); SS=length(S(1,:));
    
    %%% Pad matrices to ensure all matrices same size
    if SR<SS
        if R(1,end)==1
            C=ones(N,SS-SR);
        else
            C=zeros(N,SS-SR);
        end
        R=[R,C];
    else
        C=ones(N,SR-SS)*S(1,end);
        S=[S,C];
    end
    S=(S+R);
    I  % Display current step
    %%%%___________________________________________________________________
end
S=S./I; % Divide aggregated matrix by number of simulations to calculate state probability over time
output = sum(S(:,end))/N; %%% Calculates the fixation probability
save( 'output_bethedging.mat' , 'output' );

end


%%%% Sub functions_________________________________________________________
function [ G ] = make_random_vector( N,n )
%%% This function generates a random vector of length N with n ones and N-1
%%% zeros. If n is not an integer, generates either ceil(n) or floor(n)
%%% ones.
G=rand(N,1)<n/N;
while sum(G)~=floor(n)
    if sum(G)==ceil(n)
        break
    else
    G=rand(N,1)<n/N;
    end
end
end

function [ R ] = SimulationDBBV( T,Ap,Bp,a,c,varb,vara )
%%% This function simulates the evolution graph theory model under
%%% Death-Birth with selection on Birth dynamics. The mean fitness payoffs
%%% are specified. From this, the mean fitness of each individual is
%%% calculated, and than fitness values are randomly sampled during each
%%% time step of the simulation. 

%%% Generates an output matrix R, which captures the number of type "A"
%%% individuals at each time step of the simulation. 

%%% Generate initial matrices______________________________________________
ii=1;
R(:,ii)=Ap;
K=size(T);
N=length(T);
%%%________________________________________________________________________

%%% var_type can be changed to "Environmental" if we want to consider 
%%% between rather than within demographic variation
var_type = "Demographic"; 


while (sum(Ap)~=0 && sum(Bp)~=0)

%%%Sample fitness values___________________________________________________
    if var_type == "Environmental"
    %%Environmental Variance - individuals of a given type sample the same
    %%fitnes. Randomly sample fitness values and then build vector of
    %%individual level fitnesses.
        if vara ~= 0
            av=gamrnd(a/vara,vara);
        else
            av=a;
        end
        if varb ~= 0
            cv=gamrnd(c/varb,varb);
        else
            cv=c;
        end  
%         bv=av;
%         dv=cv;
        fA=av.*Ap;
        fB=cv.*Bp;
    else
    %%%Demographic Variance  - individuals of a given type can sample
    %%%different fitness values. Mean fitness for each individual are
    %%%calculated first, generating the vector of mean fitnesses. Values are then
    %%%randomly sampled for each index of the vector.
        fA=a.*Ap;
        fB=c.*Bp;
        if vara ~= 0
            fA=gamrnd(fA/vara,vara);
        else
            fA=fA;
        end
        if varb ~= 0
            fB=gamrnd(fB/varb,varb);
        else
            fB=fB;
        end    
        fA=fA.*Ap;
        fB=fB.*Bp;   
    end
%%%________________________________________________________________________

%%% Simulate dynamics______________________________________________________
    mm=rand; % Random number for selecting death event
    
    %%% Indexing______
    jj=1;
    %%%_______________
    while jj<=N
        Z=T(:,jj); % Subgraph (Z) of individuals connected to node jj
        if mm<=jj/N % Is node j randomly selected for death?
            fA=fA.*T(:,jj); fB=fB.*T(:,jj); % Setting fitness to zero for individuals not in the subgraph Z
            F=sum(fA(1:K(1)))+sum(fB(1:K(1))); % Total fitness across Z
            P=(fA+fB)./F; % Generate probability vector for relative fitnesses across Z
            
            %%% Indexing______
            kk=1;
            j=1;
            L=0;
            %%%_______________
            while kk<=K(2) % Generate position vector of individuals in Z
                if Z(kk)==1
                    L(j)=kk;
                    j=j+1; kk=kk+1;
                else
                    kk=kk+1;
                end
            end
            nn=rand; % Random value for selecting birth event
            
            %%% Indexing______
            ll=1;
            j=1;
            %%%_______________
            while ll<=N
                if nn<=sum(P(1:L(ll))) %Find position of individual selecting for birth
                    Ap(jj)=Ap(L(j)); % Update type of node jj to that of the selected node
                    ll=N+1; % Exit
                else
                    ll=ll+1; % Loop through position vector
                    j=j+1;
                end 
            end
            jj=K(1)+1; % Exit
        else
            jj=jj+1; % Loop over nodes until we find the one selected for death
        end

    end
    
    %%% Update vector of types B individuals____
    for i=1:K(1)
        if Ap(i)==1
            Bp(i)=0;
        else
            Bp(i)=1;
        end
    end
    %%%_________________________________________
    
    %%% Update loop index and update output matrix (R) with current system state______
    ii=ii+1;
    R(:,ii)=Ap; 
    %%%_______________________________________________________________________________
end
%%%________________________________________________________________________
%%%________________________________________________________________________

end

function [ M ] = make_complete_graph(N)

M=ones(N,N);
M=M-eye(N);

end

function [G] = make_constant_n_symmetric_matrix(N,n)
G=zeros(N);
G=sparse(G);
R=1:N;
its=0
de=1;
% while  find(sum(G^2,2)==1)~=0 

while de>0 
    its=its+1
    G_old=G;
    for i=1:N
        k=  find(G(i,:)==1);
        
        if length(k)>=n
            p=randperm(length(k));
            p=p(1:n);
            r=k(p);
            
            G(i,:)=0;
            G(:,i)=0;
            G(i,r)=1;
            G(r,i)=1;
            
        else
            A=find(G(i,:)==1);
            q=length(A);
            A=union(A,i); %cant include itself 
            
            V=setdiff(R,A);
            
            S=randperm(length(V));
            
            r=n-q; %number of entries remaining to be filled
            
            s=S(1:r);
            
            G(i,V(s))=1;
            G(V(s),i)=1;
            
        end
    end
    de=0;
%     de=abs(sum(sum(G-G_old)));
    for i=1:N
        if sum(G(i,:))~=n
            de=1;
        end
    end
     if find(sum(G^2,2)==1)~=0 
         de=1;
     else
     end
     
% end    
end
end

function [ G ] = make_scalefree_graph( N,m0,m )
%%% Generates a random scale free graph
G=make_complete_graph(m0);


for j=m0+1:N
    S=length(G);
    G(j,:)=0;
    G(:,j)=0;
    K=sum(G,2);
    L=0;
    for i=1:m
        while i>0
        t(i)=randi(S);
        l=t(i);
        if find(L==l)>0
        else
            
            L(i)=l;
        P=K(l)/sum(K);
        r=rand;
        if r<P
            G(j,l)=1;
            G(l,j)=1;
            break
        else 
        end
        end
        end
    end
end

    
    
        
       
end

function [ G ] = make_random_graph_strict( N,n )
%%% Generates an Erdős–Rényi such that the average degree is exactly equal to n

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

function [ G ] = make_random_graph( N,n )
%%% Generates an Erdős–Rényi such that the average degree is approximately equal to n

G=rand(N)<n/N;
G=tril(G);
G=G+G';
G=G-diag(diag(G));
while  find(sum(G)==0)~=0 | find(sum(G^2,2)==1)~=0 
    G=rand(N)<n/N;
    G=tril(G);
    G=G+G';
    G=G-diag(diag(G));
end
end