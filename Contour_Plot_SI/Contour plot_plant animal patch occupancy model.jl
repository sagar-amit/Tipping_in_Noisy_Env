
#=#############-----Code for generating the contour plot to investigate combined influence of q,k on the 
#tipping threshold of the control parameter in PATCH OCCUPANCY MODEL OF MUTUALISM (figure.C2(f) in SI)-----=#######

using Distributed  #package for multiprocessing in julia
nprocs()
addprocs(11) # add precessors
w=workers() #set workers


using Plots
@everywhere using Statistics
#=function to generate mean tipping threshold value across the replicates for each temporal
correlation value k and species response correlation value q=#
@everywhere function Mn(ss,k,q)
    T=150000               #Total time
    dt=0.1 #step size
    N=convert(Int,T/dt) #Number of iteration steps
    sig_X=0.03  #noise intensity for plant population
    sig_Y=0.03  ##noise intensity for animal population
    
    c_time=zeros(ss) # zero vector generated to store the tipping thresholds
    #generation of tipping threshold for each replicate
    for ll in 1:ss

        noise=zeros(N,2)
        eta = abs.(randn(2))
        for i in 1:N
            noise[i,:]=copy(eta)
            eta=k.*eta.+sqrt(1-k^2).*(sqrt(abs(q)).*randn(1).+sqrt(1-abs(q)).*randn(2))
        end 

        rept=zeros(N,2)

        rept[:,1]=sort(randn(N))
        rept[:,2]=sort(randn(N))

        n1=invperm(sortperm(noise[:,1]))
        n2=invperm(sortperm(noise[:,2]))
        noise[:,1]=rept[n1,1]         
        noise[:,2]=rept[n2,2]


        m0=0.2 # initial value of the control parameter
        x0=0.72 #plant density at T=0
        y0=0.47 #animal density at T=0



        x=zeros(N)
        y=zeros(N)
        m=zeros(N)



        x[1]=x0
        y[1]=y0
        m[1]=m0


        noise.=sqrt(dt).*noise


        for i in 2:N
            x[i]=x[i-1]+(2.0*y[i-1]*(1-m[i-1]-x[i-1])-0.1*x[i-1])*dt+x[i-1]*sig_X*noise[i-1,1]
            y[i]=y[i-1]+(2.0*y[i-1]*(x[i-1]-y[i-1])-0.5*y[i-1])*dt+y[i-1]*sig_Y*noise[i-1,2]
            m[i]=m[i-1]+(0.3/T)*dt #Control parameter gradually increases with time
        end 


        cr_v=0.18 # critical value of tipping
        for (idx,value) in enumerate(reverse(x))
            if value >= cr_v
                c_time[ll]=N-idx+1
                break
            end
        end
        


    end

    ## return mean tipping thresold value
    return 0.2.+mean(c_time).*(0.3/T)*dt
    
end


q=[0:0.1:1;]

@everywhere function task(q)
    
    k=[0:0.1:0.9;]
    sub_mat=zeros(length(k),1)
    for j in 1:length(k)
        sub_mat[j,1]=Mn(1000,k[j],q) #ss=1000 replicates generated
        println([j,q])
    end

    return sub_mat
end

#Distributed task across cores

t1=@spawnat w[1] task(q[1])
t2=@spawnat w[2] task(q[2])
t3=@spawnat w[3] task(q[3])
t4=@spawnat w[4] task(q[4])
t5=@spawnat w[5] task(q[5])
t6=@spawnat w[6] task(q[6])
t7=@spawnat w[7] task(q[7])
t8=@spawnat w[8] task(q[8])
t9=@spawnat w[9] task(q[9])
t10=@spawnat w[10] task(q[10])
t11=@spawnat w[11] task(q[11])


s1=fetch(t1)
s2=fetch(t2)
s3=fetch(t3)
s4=fetch(t4)
s5=fetch(t5)
s6=fetch(t6)
s7=fetch(t7)
s8=fetch(t8)
s9=fetch(t9)
s10=fetch(t10)
s11=fetch(t11)


mat=hcat(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11)

## contour plot
contourf([0:0.1:0.9;],[0:0.1:1;],transpose(mat),lw=0.1,levels=15,color=:turbo,xlab="Temporal correlation",
ylab="Species response correlation") 

