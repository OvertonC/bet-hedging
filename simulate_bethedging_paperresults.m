%%%% Main function_________________________________________________________
function simulate_bethedging_paperresults
load input1.mat out; %%% input file contains matrices for degree variation 
                     %%% and seed for random number generator
                     %%% by generating multiple input files with different 
                     %%% seeds we can run the model in parallel
                     %%% and aggregate the results at the end.
input = out; 
rng(input.seed); % set random seed


%%% Define mean fitness values
a=0.95; % Mean bet-hedger fitness
c=1; % Mean normal-type fitness


%%% Specify number of simulations to be run. Better to use low number and
%%% run in parallel.
%%% To generate the "Degree variation" results in the paper, we used
%%% 30,000,000 simulation, run in parallel. We set tot = 300 and ran in 
%%% parallel using 100000 different seeds
%%% To generate the "Size" and "Degree variation" results in the paper, 
%%% we used 1,000,000 simulations. For "Size", we set tot = 10 and run in
%%% parallel using 100000 different seeds. For "Average degree" we set tot
%%% = 1000 and run in parallel using 1000 different seeds.
%%% The procedure was setup in this way to keep the run time of each batch
%%% as low as possible (under 20 minutes) whilst using as few batches as
%%% possible. This work was facilitated through the use of the ARC Condor
%%% High Throughput Computing system at the University of Liverool, which
%%% allowed 3000 jobs to be run simultaneously. 

tot = 10;


%%% What result do we want to generate?
result_type = "Average degree"; % "Degree variation"; "Size"; "Average degree";

if result_type == "Degree variation"
    index_limit = 16;
    for k=1:index_limit
        inname = ['matrix',num2str(k)];
        T = input.(inname);    
        [results,results1,results2] = loop_variation1(tot,T,a,c);
        outputresults_degreevariation(k,:) = [results,results1,results2];
    end
    save( 'output_degreevariation_test.mat' , 'outputresults_degreevariation' );
elseif result_type == "Size"
    index_limit = 3;
    for k=1:index_limit
        if k==1
            T=make_constant_n_symmetric_matrix(50,4);
        elseif k==2
            T=make_constant_n_symmetric_matrix(100,4);
        elseif k==3
            T=make_constant_n_symmetric_matrix(200,4);
        else
            print("Population size not set up")
        end
        [results,results1,results2] = loop_variation2(tot,T,a,c);
        outputresults_popsize(k,:) = [results,results1,results2];
    end
    save( 'output_popsize_test.mat' , 'outputresults_popsize' );   
elseif result_type == "Average degree"
    index_limit = 3;
    for k=1:index_limit
        if k==1
            T=make_constant_n_symmetric_matrix(50,4);
        elseif k==2
            T=make_constant_n_symmetric_matrix(50,8);
        elseif k==3
            T=make_constant_n_symmetric_matrix(50,16);
        elseif k ==4
            T=make_complete_graph(50);
        else
            print("Average degree not set up")
        end
        for z = 1:3
            if z == 1
                [results,results1,results2] = loop_variation3(tot,T,a,c);
                outputresults_average_0(k,:) = [results,results1,results2];
            elseif z == 2
                [results,results1,results2] = loop_variation4(tot,T,a,c); 
                outputresults_average_1(k,:) = [results,results1,results2];
            else
                [results,results1,results2] = loop_variation5(tot,T,a,c);
                outputresults_average_2(k,:) = [results,results1,results2];
            end
        end        
    end
    save( 'output_average_0_test.mat' , 'outputresults_average_0' );
    save( 'output_average_1_test.mat' , 'outputresults_average_1' );
    save( 'output_average_2_test.mat' , 'outputresults_average_2' );
end
end


%%%% Sub functions_________________________________________________________
function [results,results1,results2] = loop_variation1(tot,T,a,c)
for g = 1:10
    vara = 2.1053;
    varb = 2.2 + 0.04*(g-1);       
    [R1,R2,R3] = call_simulator(tot,T,a,c,varb,vara);
    results(g) = R3;
    results1(g) = R1;
    results2(g) = R2;
end
end

function [results,results1,results2] = loop_variation2(tot,T,a,c)
for g = 1:7
    vara = 2.1053;
    varb = 2.1 + 0.1*(g-1);       
    [R1,R2,R3] = call_simulator(tot,T,a,c,varb,vara);
    results(g) = R3;
    results1(g) = R1;
    results2(g) = R2;
end
end

function [results,results1,results2] = loop_variation3(tot,T,a,c)
for g = 1:11
    vara = 0;
    varb = 0.+0.1*(g-1);       
    [R1,R2,R3] = call_simulator(tot,T,a,c,varb,vara);
    results(g) = R3;
    results1(g) = R1;
    results2(g) = R2;
end
end

function [results,results1,results2] = loop_variation4(tot,T,a,c)
for g = 1:11
    vara = 1.0526;
    varb = 1+0.1*(g-1);       
    [R1,R2,R3] = call_simulator(tot,T,a,c,varb,vara);
    results(g) = R3;
    results1(g) = R1;
    results2(g) = R2;
end
end

function [results,results1,results2] = loop_variation5(tot,T,a,c)
for g = 1:11
    vara = 2.1053;
    varb = 2.1+0.01*(g-1);       
    [R1,R2,R3] = call_simulator(tot,T,a,c,varb,vara);
    results(g) = R3;
    results1(g) = R1;
    results2(g) = R2;
end
end

function [ G ] = make_random_vector2( N,n )
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

function [R1, R2, R3] = call_simulator(tot,T,a,c,varb,vara)
N = size(T,2); 
S=zeros(N,N);
R=zeros(N,N);
for indexj=1:2
    if indexj==1
        for I=1:tot
            N = size(T,2);
            Ap=make_random_vector2(N,1); Bp=1-Ap;
            [R]=SimulationDBBV(T,Ap,Bp,a,c,varb,vara);
            SR=length(R(1,:)); SS=length(S(1,:));
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
            I
        end
        S=S./I;
        output(indexj)=sum(S(:,end))/N;
        S=zeros(N,N);
        R=zeros(N,N);
    else
        for I=1:tot
            N = size(T,2);
            Ap=make_random_vector2(N,N-1); Bp=1-Ap;
            [R]=SimulationDBBV(T,Ap,Bp,a,c,varb,vara);
            SR=length(R(1,:)); SS=length(S(1,:));
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
            I
        end
        S=S./I;
        output(indexj)=sum(S(:,end))/N;
        S=zeros(N,N);
        R=zeros(N,N);
    end
end
R3 = output(1)/(1-output(2));
R1 = output(1);
R2 = output(2);
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


