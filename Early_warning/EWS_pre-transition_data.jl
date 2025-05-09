#################### CODE TO GENERATE DATA FOR CALCULATION OF EARLY WARNING SIGNALS ##############################

# For each model the pre-transition time series data is generated for one species (x) (species 1)
# with this code one can also generate the pre-transition time series data for the another species (y) (species 2)
##################################################################################################################

using CSV
using DataFrames
using Statistics

####################################### FOR FIGURES IN THE MAIN TEXT ##########################################################



#################################################################################################
#################################################################################################
#################### DATA GENERATION FOR SHALLOW-LAKE MODEL######################################
#################################################################################################
#################################################################################################


T=20000 #total time
dt=0.5 #step size
N=convert(Int,T/dt)  #Number of iteration steps

## Function to generate ss stochastic trajectories for each q
function Mn(q,ss)
  sig_X=0.02 #Noise intensity for Macrophyte population
  sig_Y=0.02 #Noise intensity for Algae population
  k=0.0 ## temporal correlation fixed at k=0.0

  gen_noise = [zeros(N,2) for _ in 1:ss]
  noise=zeros(N,2)
  eta = abs.(randn(2))
  
  ## Noise generation
  for j in 1:ss
    for i in 1:N
      noise[i,:]=copy(eta)
      eta=k.*eta.+sqrt(1-k^2).*(sqrt(abs(q)).*randn(1).+sqrt(1-abs(q)).*randn(2))
    end
    
    gen_noise[j]=copy(noise)
  end
   
  #Spectral mimicry
  rept=zeros(N,2)
  gen_rept= [zeros(N,2) for _ in 1:ss]
    
  for j in 1:ss
    rept[:,1]=sort(randn(N))
    rept[:,2]=sort(randn(N))
    gen_rept[j]=copy(rept)
  end
  
  for j in 1:ss
    n1=invperm(sortperm(gen_noise[j][:,1]))
    n2=invperm(sortperm(gen_noise[j][:,2]))
    gen_rept[j][:,1]=gen_rept[j][n1,1]
    gen_rept[j][:,2]=gen_rept[j][n2,2]
  end
  
  #parameter values
  r1=0.05
  r2=0.1
  K=1
  h1=0.2
  h2=2
  p=4
  E0=4 # control parameter (Nutrient loading)
  
  #Initial values
  x0=0.987
  y0=0.674

  x=zeros(N,ss)
  y=zeros(N,ss)
  E=zeros(N,ss)

  x[1,:].=x0
  y[1,:].=y0
  E[1,:].=E0
    
  for j in 1:ss
    gen_rept[j]=sqrt(dt).*gen_rept[j]
  end
  
  for j in 1:ss
    
    for i in 2:N
      x[i,j]=x[i-1,j]+r1*x[i-1,j]*(1-(x[i-1,j]/(K*h2^p))*(h2^p+y[i-1,j]^p))*dt+x[i-1,j]*sig_X*gen_rept[j][i-1,1]
      y[i,j]=y[i-1,j]+r2*y[i-1,j]*(1-(y[i-1,j]/(h1*E[i-1,j]))*(h1+x[i-1,j]))*dt+y[i-1,j]*sig_Y*gen_rept[j][i-1,2]
      E[i,j]=E[i-1,j]+(4/T)*dt
    end

  end

  return x 
  #for the another species (algae) return y
end


ss=150 ## number of replicates for each q
q=[0:0.225:0.9;] # species response correlation values: q=0.0, 0.225, 0.45, 0.675, 0.9

# store the matrices generated for each q as a list, where each matrix consists the 150 trajectories in its coulmns
# for each q and so there will a total of 750 trajectories in the list
rep_mat=[zeros(N,ss) for _ in 1:length(q)]

for i in 1:length(q)
  rep_mat[i]=Mn(q[i],ss)
end

lind=zeros(Int,ss,length(q))
## sotre the index where transition occurs
for i in 1:length(q)
  for j in 1:ss
    cr_v=0.375    ##transition value for tipping  (we have to consider cr_v=2.713 for algae (y))
    for (idx,value) in enumerate(reverse(rep_mat[i][:,j]))
      if value >= cr_v  ## set, value <= cr_v for algae (y)
        lind[j,i]=trunc.(Int,N-idx+1)
        break
      end
    end
  end
end

B=minimum(lind) # miminim of all index where transition occurs

gen_mat=[zeros(B,ss) for _ in 1:length(q)]
## store pre-transition time series  by reducing each population time series to the length 
#of the shortest time series where transition occurs across all 750 simulated dataset.
for i in 1:length(q)
  for j in 1:ss
    gen_mat[i][:,j]=rep_mat[i][1:minimum(lind),j]
  end
end

# save the generated pre-transition data sets for Macrophyte (x) which being used to calculate EWS  (in similar way one can generate the pre-transition data for Algae (y))

data=DataFrame(gen_mat[1],:auto)
CSV.write("lake150_M_ews_q0.0_k0.0.csv",data)

data=DataFrame(gen_mat[2],:auto)
CSV.write("lake150_M_ews_q0.225_k0.0.csv",data)

data=DataFrame(gen_mat[3],:auto)
CSV.write("lake150_M_ews_q0.45_k0.0.csv",data)

data=DataFrame(gen_mat[4],:auto)
CSV.write("lake150_M_ews_q0.675_k0.0.csv",data)

data=DataFrame(gen_mat[5],:auto)
CSV.write("lake150_M_ews_q0.9_k0.0.csv",data)



######################################################################################
###################### DATA GENERATION FOR PREDATOR PREY MODEL########################
######################################################################################

T=20000 ##Total time
dt=0.5  # step size
N=convert(Int,T/dt) #Number of iteration steps

function Mn(q,ss)#= function to calculate time series, where q represents the degree of species
                   response correlation and ss is the number of replicas for each q =#
  
  k=0.0 ## temporal correlation fixed at k=0.0
  sig_X=0.013 #noise intensity for prey population
  sig_Y=0.013 ##noise intensity for predator population

  gen_noise = [zeros(N,2) for _ in 1:ss] # list of zero matrices to store original pairwise noise
  noise=zeros(N,2) # zero matrix used to generate noise
  eta = abs.(randn(2)) # initial values of noise
  ###generating the pairwise noise
  for j in 1:ss
    for i in 1:N
      noise[i,:]=copy(eta)
      eta=k.*eta.+sqrt(1-k^2).*(sqrt(abs(q)).*randn(1).+sqrt(1-abs(q)).*randn(2))
    end 
    gen_noise[j]=copy(noise)
  end
  ###spectral mimicry used to generate noise from the original one so that it follows standard normal distribution  
  rept=zeros(N,2)
  gen_rept= [zeros(N,2) for _ in 1:ss]

  for j in 1:ss
    rept[:,1]=sort(randn(N))
    rept[:,2]=sort(randn(N))
    gen_rept[j]=copy(rept)
  end
  
  ####### replace the generated noise with standard normal noise 
  for j in 1:ss
    n1=invperm(sortperm(gen_noise[j][:,1]))
    n2=invperm(sortperm(gen_noise[j][:,2]))
    gen_rept[j][:,1]=gen_rept[j][n1,1]
    gen_rept[j][:,2]=gen_rept[j][n2,2]
      
  end
  #### spectral mimicry ends here 
    
  
  r=0.2 #Intrinsic growth rate of prey
  a=0.4 #attack rate
  b=7    #saturation constant
  k1=1.2 #carrying capacity for the predator populatiion
  d=6   #conversion of biomass from predator to prey (d=c*a)
  m0=0.01  #initial value of the control parameter (predator mortality rate)
   
   
  x0=0.003 #initial prey density at T=0
  y0=0.509 #initial predator density at T=0

  x=zeros(N,ss)
  y=zeros(N,ss)
  m=zeros(N,ss)

  x[1,:].=x0
  y[1,:].=y0
  m[1,:].=m0

  for j in 1:ss
    gen_rept[j]=sqrt(dt).*gen_rept[j]  
  end

  for j in 1:ss

    for i in 2:N
      x[i,j]=x[i-1,j].+(r.*x[i-1,j].*(1-x[i-1,j]).-((a.*x[i-1,j].*y[i-1,j])./(1 .+b.*x[i-1,j]))).*dt.+x[i-1,j].*sig_X.*gen_rept[j][i-1,1]
      y[i,j]=y[i-1,j].+(((d.*x[i-1,j].*y[i-1,j])./(1 .+b.*x[i-1,j])).*(1-y[i-1,j]./k1)-m[i-1,j].*y[i-1,j]).*dt.+y[i-1,j].*sig_Y.*gen_rept[j][i-1,2]
      m[i,j]=m[i-1,j].+(0.14/T)*dt
    end 
     
  end
    
  return x   #returning the time series corresponding to prey population, in case of predator population return y
end


ss=150  ## number of replicas generated for each q
q=[0:0.225:0.9;] ## species response correlation values: 0.0, 0.225, 0.45, 0.675, 0.9



rep_mat=[zeros(N,ss) for _ in 1:length(q)] # list of zero vectors to store time series data for each q
for i in 1:length(q)
  rep_mat[i]=Mn(q[i],ss)
end



lind=zeros(Int,ss,length(q)) # zero matrix to store the time index at which transition occurs where coulmns store the time index for ss replicas for each q
for i in 1:length(q)
  for j in 1:ss
    cr_v=0.265    #transition value for tipping, we have to consider cr_v=0.961 for predator population (y)
    for (idx,value) in enumerate(reverse(rep_mat[i][:,j]))
      if value <= cr_v 
        lind[j,i]=trunc.(Int,N-idx+1)
        break
      end
    end
  end
end


B=minimum(lind)
gen_mat=[zeros(B,ss) for _ in 1:length(q)] #list of zero vectors to store pre-transition time series data
                                      
for i in 1:length(q)
  for j in 1:ss 
    gen_mat[i][:,j]=rep_mat[i][1:minimum(lind),j]
  end
end

# save the generated data sets for x which being used to calculate EWS and in similar one can generate the 
#pre-transition data for y

data=DataFrame(gen_mat[1],:auto)
CSV.write("pp150_x_ews_q0.0_k0.0.csv",data) # 150 replicas for q=0.0

data=DataFrame(gen_mat[2],:auto)
CSV.write("pp150_x_ews_q0.225_k0.0.csv",data) # 150 replicas for q=0.225

data=DataFrame(gen_mat[3],:auto)
CSV.write("pp150_x_ews_q0.45_k0.0.csv",data) # 150 replicas for q=0.45

data=DataFrame(gen_mat[4],:auto)
CSV.write("pp150_x_ews_q0.675_k0.0.csv",data) # 150 replicas for q=0.675

data=DataFrame(gen_mat[5],:auto)
CSV.write("pp150_x_ews_q0.9_k0.0.csv",data) # 150 replicas for q=0.9



####################################################################################
####################################################################################
############# DATA GENERATION FOR GENERAL MODEL OF MUTUALISM #######################
####################################################################################
####################################################################################

T=4000
dt=0.1
N=convert(Int,T/dt)

function Mn(q,ss)
  sig_X=0.03
  sig_Y=0.03
  k=0.0
  
  gen_noise = [zeros(N,2) for _ in 1:ss]
  noise=zeros(N,2)
  eta = abs.(randn(2))
    
  for j in 1:ss
    for i in 1:N
      noise[i,:]=copy(eta)
      eta=k.*eta.+sqrt(1-k^2).*(sqrt(abs(q)).*randn(1).+sqrt(1-abs(q)).*randn(2))
    end
    gen_noise[j]=copy(noise)
  end
    
  rept=zeros(N,2)
  gen_rept= [zeros(N,2) for _ in 1:ss]
    
  for j in 1:ss
    rept[:,1]=sort(randn(N))
    rept[:,2]=sort(randn(N))
    gen_rept[j]=copy(rept)
  end
    
  for j in 1:ss
    n1=invperm(sortperm(gen_noise[j][:,1]))
    n2=invperm(sortperm(gen_noise[j][:,2]))
    gen_rept[j][:,1]=gen_rept[j][n1,1]
    gen_rept[j][:,2]=gen_rept[j][n2,2]
  end
    
     
  r1=-1
  r2=-0.25
  a1=2
  a2=2
  b1=2
  b2=2
  c1=1
  c2=1
  h1=0.5
  h2=0.5

  x0=1.79
  y0=2.318
  
  x=zeros(N,ss)
  y=zeros(N,ss)
  a=zeros(N,ss)

  x[1,:].=x0
  y[1,:].=y0
  a[1,:].=r2
  
  for j in 1:ss
    gen_rept[j]=sqrt(dt).*gen_rept[j]
  end
  
  for j in 1:ss
    for i in 2:N
      x[i,j]=x[i-1,j].+(x[i-1,j].*(r1.-c1.*x[i-1,j].+b1.*(a1.*y[i-1,j]./(1.0.+a1.*h1.*y[i-1,j])))).*dt.+x[i-1,j].*sig_X.*gen_rept[j][i-1,1]
      y[i,j]=y[i-1,j].+(y[i-1,j].*(a[i-1,j].-c2.*y[i-1,j].+b2.*(a2.*x[i-1,j]./(1.0.+a2.*h2.*x[i-1,j])))).*dt.+y[i-1,j].*sig_Y.*gen_rept[j][i-1,2]
      a[i,j]=a[i-1,j].-((1.15-0.25)/T)*dt
    end
  end
    
  return x
end


ss=150
q=[0:0.225:0.9;]

rep_mat=[zeros(N,ss) for _ in 1:length(q)]
for i in 1:length(q)
  rep_mat[i]=Mn(q[i],ss)
end


lind=zeros(Int,ss,length(q))

for i in 1:length(q)
  for j in 1:ss
    cr_v=0.5  # transition value is also the same for mutualist 2 (y)
    for (idx,value) in enumerate(reverse(rep_mat[i][:,j]))
      if value >= cr_v #same for mutualist 2 (y)
        lind[j,i]=trunc.(Int,N-idx+1)
        break
      end
    end
  end
end

B=minimum(lind)

gen_mat=[zeros(B,ss) for _ in 1:length(q)]

for i in 1:length(q)
  for j in 1:ss
    gen_mat[i][:,j]=rep_mat[i][1:minimum(lind),j]
  end
end

## save data sets being analyzed for EWSs corresponding to mutualist 1 population (x)
data=DataFrame(gen_mat[1],:auto)
CSV.write("mutob150_M1_ews_q0.0_k0.0.csv",data)

data=DataFrame(gen_mat[2],:auto)
CSV.write("mutob150_M1_ews_q0.225_k0.0.csv",data)

data=DataFrame(gen_mat[3],:auto)
CSV.write("mutob150_M1_ews_q0.45_k0.0.csv",data)

data=DataFrame(gen_mat[4],:auto)
CSV.write("mutob150_M1_ews_q0.675_k0.0.csv",data)

data=DataFrame(gen_mat[5],:auto)
CSV.write("mutob150_M1_ews_q0.9_k0.0.csv",data)



####################################### FOR FIGURES IN THE SUPPLEMENTERY INFORMATION ##########################################################


##################################################################################
##################################################################################
############### DATA GENERATION FOR MACROALGAE CORAL MODEL########################
##################################################################################
##################################################################################

T=20000
dt=0.5
N=convert(Int,T/dt)


function Mn(q,ss)
  sig_X=0.02
  sig_Y=0.02
  k=0.0
  
  
  gen_noise = [zeros(N,2) for _ in 1:ss]
  noise=zeros(N,2)
  eta = abs.(randn(2))
    
  for j in 1:ss
    for i in 1:N
      noise[i,:]=copy(eta)
      eta=k.*eta.+sqrt(1-k^2).*(sqrt(abs(q)).*randn(1).+sqrt(1-abs(q)).*randn(2))
    end
        
    gen_noise[j]=copy(noise)
  end
    
  rept=zeros(N,2)
  gen_rept= [zeros(N,2) for _ in 1:ss]
    
  for j in 1:ss
    rept[:,1]=sort(randn(N))
    rept[:,2]=sort(randn(N))
    gen_rept[j]=copy(rept)
  end
    
  for j in 1:ss
    n1=invperm(sortperm(gen_noise[j][:,1]))
    n2=invperm(sortperm(gen_noise[j][:,2]))
    gen_rept[j][:,1]=gen_rept[j][n1,1]
    gen_rept[j][:,2]=gen_rept[j][n2,2]
  end
     
  m0=0.53 
  x0=0.0025
  y0=0.5606

  x=zeros(N,ss)
  y=zeros(N,ss)
  m=zeros(N,ss)

  x[1,:].=x0
  y[1,:].=y0
  m[1,:].=m0
    
  for j in 1:ss
    gen_rept[j]=sqrt(dt).*gen_rept[j]  
  end

  for j in 1:ss
    for i in 2:N
      x[i,j]=x[i-1,j].+(x[i-1,j].*(0.1.*y[i-1,j].-(m[i-1,j]./(1-y[i-1,j])).+0.77.*(1.0.-x[i-1,j].-y[i-1,j])).+0.005.*(1-x[i-1,j]-y[i-1,j])).*dt+x[i-1,j].*sig_X.*gen_rept[j][i-1,1]
      y[i,j]=y[i-1,j].+(y[i-1,j].*(0.55.*(1.0.-x[i-1,j].-y[i-1,j]).-0.24.-0.1.*x[i-1,j])).*dt.+y[i-1,j].*sig_Y.*gen_rept[j][i-1,2]
      m[i,j]=m[i-1,j]-((0.53-0.15)/T)*dt
    end  
  end
    
  return x ## return y for coral

end


ss=150
q=[0:0.225:0.9;]

rep_mat=[zeros(N,ss) for _ in 1:length(q)]

for i in 1:length(q)
  rep_mat[i]=Mn(q[i],ss)
end



lind=zeros(Int,ss,length(q))

for i in 1:length(q)
  for j in 1:ss
    cr_v=0.26 ##we have to consider cr_v=0.255 for coral (y)
    for (idx,value) in enumerate(reverse(rep_mat[i][:,j]))
      if value <= cr_v # set value >= cr_v for coral (y)
        lind[j,i]=trunc.(Int,N-idx+1)
        break
       end
    end
  end
end

B=minimum(lind)

gen_mat=[zeros(B,ss) for _ in 1:length(q)]

for i in 1:length(q)
  for j in 1:ss
    gen_mat[i][:,j]=rep_mat[i][1:minimum(lind),j]
  end
end

## save data sets being analyzed for EWSs corresponding to MACROALGAE population
data=DataFrame(gen_mat[1],:auto)
CSV.write("Maccoral150_M_ews_q0.0_k0.0.csv",data)

data=DataFrame(gen_mat[2],:auto)
CSV.write("Maccoral150_M_ews_q0.225_k0.0.csv",data)

data=DataFrame(gen_mat[3],:auto)
CSV.write("Maccoral150_M_ews_q0.45_k0.0.csv",data)

data=DataFrame(gen_mat[4],:auto)
CSV.write("Maccoral150_M_ews_q0.675_k0.0.csv",data)

data=DataFrame(gen_mat[5],:auto)
CSV.write("Maccoral150_M_ews_q0.9_k0.0.csv",data)



###########################################################################################
###########################################################################################
############ DATA GENERATION FOR PLANT-POLLINATOR MODEL OF MUTUALISM########################
###########################################################################################
###########################################################################################

T=20000
dt=0.5
N=convert(Int,T/dt)

function Mn(q,ss)
  
  k=0.0
  sig_X=0.006
  sig_Y=0.006
  
  gen_noise = [zeros(N,2) for _ in 1:ss]
  noise=zeros(N,2)
  eta = abs.(randn(2))
  for j in 1:ss
    for i in 1:N
      noise[i,:]=copy(eta)
      eta=k.*eta.+sqrt(1-k^2).*(sqrt(abs(q)).*randn(1).+sqrt(1-abs(q)).*randn(2))
    end 
    gen_noise[j]=copy(noise)
  end
    
  rept=zeros(N,2)
  gen_rept= [zeros(N,2) for _ in 1:ss]
  for j in 1:ss
    rept[:,1]=sort(randn(N))
    rept[:,2]=sort(randn(N))
    gen_rept[j]=copy(rept)
        
  end
  for j in 1:ss
    n1=invperm(sortperm(gen_noise[j][:,1]))
    n2=invperm(sortperm(gen_noise[j][:,2]))
    gen_rept[j][:,1]=gen_rept[j][n1,1]
    gen_rept[j][:,2]=gen_rept[j][n2,2]
      
  end
    
 a0=1.8 
 x0=0.108
 y0=0.2
 
 x=zeros(N,ss)
 y=zeros(N,ss)
 a=zeros(N,ss)

 x[1,:].=x0
 y[1,:].=y0
 a[1,:].=a0

  for j in 1:ss
   gen_rept[j]=sqrt(dt).*gen_rept[j]
  end

 for j in 1:ss
   for i in 2:N
     x[i,j]=x[i-1,j].+(x[i-1,j].*(1.42.*(0.5.+0.5.*(a[i-1,j].*x[i-1,j].*y[i-1,j]./(1.0.+a[i-1,j].*0.01.*x[i-1,j].+a[i-1,j].*y[i-1,j].*x[i-1,j]))).-0.61.*x[i-1,j].-0.67)).*dt.+x[i-1,j].*sig_X.*gen_rept[j][i-1,1]
     y[i,j]=y[i-1,j].+(y[i-1,j].*(1.0.+a[i-1,j].*x[i-1,j]./(1.0.+a[i-1,j].*0.01.*x[i-1,j]).-2.0.*y[i-1,j].-0.8)).*dt.+y[i-1,j].*sig_Y.*gen_rept[j][i-1,2]
     a[i,j]=a[i-1,j].+((2.12-1.8)/T)*dt
    end
  end 
  return x  #return y for animal population
end


ss=150
q=[0:0.225:0.9;]
rep_mat=[zeros(N,ss) for _ in 1:length(q)]
for i in 1:length(q)
  rep_mat[i]=Mn(q[i],ss)
end


lind=zeros(Int,ss,length(q))
for i in 1:length(q)
  for j in 1:ss
    cr_v=0.345 ##transition value; set 0.45 for animal population
    for (idx,value) in enumerate(reverse(rep_mat[i][:,j]))
      if value <= cr_v
        lind[j,i]=trunc.(Int,N-idx+1)
        break
      end
    end
  end
end

B=minimum(lind)

gen_mat=[zeros(B,ss) for _ in 1:length(q)]

for i in 1:length(q)
  for j in 1:ss 
    gen_mat[i][:,j]=rep_mat[i][1:minimum(lind),j]
  end
end

## save data sets being analyzed for EWSs corresponding to plant population (x)
data=DataFrame(gen_mat[1],:auto)
CSV.write("mutbfac150_P_ews_q0.0_k0.0.csv",data)

data=DataFrame(gen_mat[2],:auto)
CSV.write("mutbfac150_P_ews_q0.225_k0.0.csv",data)

data=DataFrame(gen_mat[3],:auto)
CSV.write("mutbfac150_P_ews_q0.45_k0.0.csv",data)

data=DataFrame(gen_mat[4],:auto)
CSV.write("mutbfac150_P_ews_q0.675_k0.0.csv",data)

data=DataFrame(gen_mat[5],:auto)
CSV.write("mutbfac150_P_ews_q0.9_k0.0.csv",data)

########################################################################################
########################################################################################
############### DATA GENERATION FOR PATCH OCCUPANCY MODEL OF MUTUALISM###################
#########################################################################################
#########################################################################################

T=20000
dt=0.5
N=convert(Int,T/dt)

function Mn(q,ss)
 
  sig_X=0.03
  sig_Y=0.03
  k=0.0
  
  gen_noise = [zeros(N,2) for _ in 1:ss]
  noise=zeros(N,2)
  eta = abs.(randn(2))
  for j in 1:ss
    for i in 1:N
      noise[i,:]=copy(eta)
      eta=k.*eta.+sqrt(1-k^2).*(sqrt(abs(q)).*randn(1).+sqrt(1-abs(q)).*randn(2))
    end 
    gen_noise[j]=copy(noise)
  end
    
  rept=zeros(N,2)
  gen_rept= [zeros(N,2) for _ in 1:ss]
  for j in 1:ss
    rept[:,1]=sort(randn(N))
    rept[:,2]=sort(randn(N))
    gen_rept[j]=copy(rept)
        
  end
  for j in 1:ss
    n1=invperm(sortperm(gen_noise[j][:,1]))
    n2=invperm(sortperm(gen_noise[j][:,2]))
    gen_rept[j][:,1]=gen_rept[j][n1,1]
    gen_rept[j][:,2]=gen_rept[j][n2,2]
      
  end
      
  
  m0=0.2
  x0=0.723
  y0=0.473

  x=zeros(N,ss)
  y=zeros(N,ss)
  m=zeros(N,ss)

  x[1,:].=x0
  y[1,:].=y0
  m[1,:].=m0

  for j in 1:ss
    gen_rept[j]=sqrt(dt).*gen_rept[j]  
  end

  for j in 1:ss
    for i in 2:N
      x[i,j]=x[i-1,j]+(2.0*y[i-1,j]*(1-m[i-1,j]-x[i-1,j])-0.1*x[i-1,j])*dt+x[i-1,j]*sig_X*gen_rept[j][i-1,1];
      y[i,j]=y[i-1,j]+(2.0*y[i-1,j]*(x[i-1,j]-y[i-1,j])-0.5*y[i-1,j])*dt+y[i-1,j]*sig_Y*gen_rept[j][i-1,2];
      m[i,j]=m[i-1,j]+(0.3/T)*dt
    end 
  end
    
  return x # return y for fraction of patch occupied by animals
end


ss=150
q=[0:0.225:0.9;]

rep_mat=[zeros(N,ss) for _ in 1:length(q)]

for i in 1:length(q)
  rep_mat[i]=Mn(q[i],ss)
end


lind=zeros(Int,ss,length(q))
for i in 1:length(q)
  for j in 1:ss
    cr_v=0.18 ## set cr_v=0.056 for fraction of patch occupied by animals (y)
    for (idx,value) in enumerate(reverse(rep_mat[i][:,j]))
      if value >= cr_v # same for y
        lind[j,i]=trunc.(Int,N-idx+1)
        break
      end
    end
  end
end

B=minimum(lind)
gen_mat=[zeros(B,ss) for _ in 1:length(q)]

for i in 1:length(q)
  for j in 1:ss 
    gen_mat[i][:,j]=rep_mat[i][1:minimum(lind),j]
  end
end

## save data sets being analyzed for EWSs corresponding to plants occupied patches (x)
data=DataFrame(gen_mat[1],:auto)
CSV.write("mutsp150_p_ews_q0.0_k0.0.csv",data)

data=DataFrame(gen_mat[2],:auto)
CSV.write("mutsp150_p_ews_q0.225_k0.0.csv",data)

data=DataFrame(gen_mat[3],:auto)
CSV.write("mutsp150_p_ews_q0.45_k0.0.csv",data)

data=DataFrame(gen_mat[4],:auto)
CSV.write("mutsp150_p_ews_q0.675_k0.0.csv",data)

data=DataFrame(gen_mat[5],:auto)
CSV.write("mutsp150_p_ews_q0.9_k0.0.csv",data)
