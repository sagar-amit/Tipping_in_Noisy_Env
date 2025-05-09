################## Contour plot for shallow-lake model (figure.4(a)) ###############################
#################################################################################

using Distributed  #package for multiprocessing in julia
nprocs()
addprocs(11) # add precessors
w=workers() #set workers



using Plots
@everywhere using Statistics
#=function to generate mean tipping threshold value across the ss replicates for each temporal
correlation value k and species response correlation value q=#
@everywhere function Mn(ss,k,q)
  T=150000  #Total time
  dt=0.1 #step size
  N=convert(Int,T/dt) #Number of iteration steps
  sig_X=0.02 #noise intensity for Macrophyte population
  sig_Y=0.02  ##noise intensity for Algae population
    
  c_time=zeros(ss) # zero vector generated to store the tipping thresholds

  #generation of tipping threshold for each replicate
  for ll in 1:ss
    # Noise generation with spectral mimicry 
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
    
    # Parameter values
    r1=0.05
    r2=0.1
    K=1
    h1=0.2
    h2=2
    p=4
    E0=4 #initial value of the control parameter (Nutrient loading)
    
    
    # initial state
    x0=0.987
    y0=0.674
     

    x=zeros(N)
    y=zeros(N)
    E=zeros(N)



    x[1]=x0
    y[1]=y0
    E[1]=E0


    noise.=sqrt(dt).*noise

    for i in 2:N
      x[i]=x[i-1].+r1 .*x[i-1].*(1 .-(x[i-1]./(K.*h2.^p)).*(h2.^p+y[i-1].^p)).*dt.+x[i-1].*sig_X.*noise[i-1,1]
      y[i]=y[i-1].+r2.*y[i-1].*(1 .-(y[i-1]./(h1 .*E[i-1])).*(h1 .+x[i-1])).*dt.+y[i-1].*sig_Y.*noise[i-1,2]
      E[i]=E[i-1].+((4)./T).*dt 
    end 

    # Determination of tipping threshold for each replicas

    cr_v=0.375 # critical value of tipping
    for (idx,value) in enumerate(reverse(x))
      if value >= cr_v
        c_time[ll]=N-idx+1
        break
      end
    end

  end

  ## return mean tipping thresold value
  return 4.0.+mean(c_time).*((4)/T)*dt

end



q=[0:0.1:1;] # species response correlation values

# generate a function to run the simulation for each q parallely in 11 cores
@everywhere function task(q)

  k=[0:0.1:0.9;]  #temporal correlation values
  sub_mat=zeros(length(k),1)
  
  for j in 1:length(k)
    sub_mat[j,1]=Mn(1000,k[j],q) #mean tipping of ss=1000 replicates generated
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


contourf([0:0.1:0.9;],[0:0.1:1;],transpose(mat),lw=0.1,levels=15,color=:turbo, xlab="Temporal correlation",
ylab="Species response correlation") ## contour plot


