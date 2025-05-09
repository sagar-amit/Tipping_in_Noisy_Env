using DelimitedFiles
using Statistics

#q degree of species response correlation
#k degree of temporal environmental correlation
#ss number of replicate

#############################################################################################
#############################################################################################
############  GENERATION OF FLUCTUATION SIZE FOR PREDATOR-PREY MODEL ########################
#############################################################################################
#############################################################################################

function fluc(k,q,ss)
  T=1500000  # Total time
  dt=0.1     # step size
  N=convert(Int,T/dt) #Numer of itterative steps
  
  sig_X=0.013 #noise intensity in prey growth
  sig_Y=0.013 #noise intensity in predator growth
  gen_noise = [zeros(N,2) for _ in 1:ss] 
  noise=zeros(N,2)
  eta = abs.(randn(2))
  
  for j in 1:ss
    for i in 1:N 
      noise[i,:]=copy(eta)
      eta=k.*eta.+sqrt(1-k^2).*(sqrt(abs(q)).*randn(1).+sqrt(1-abs(q)).*randn(2)) #noise generating function
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

 ########### replace the generated noise with standard normal noise  
  for j in 1:ss
    n1=invperm(sortperm(gen_noise[j][:,1]))
    n2=invperm(sortperm(gen_noise[j][:,2]))
    gen_rept[j][:,1]=gen_rept[j][n1,1]
    gen_rept[j][:,2]=gen_rept[j][n2,2]
  end
  #### spectral mimicry ends here 
  
  r=0.2 #growth rate of prey
  a=0.4 #attack rate
  b=7 #saturation constant
  k1=1.2 #carrying capacity for predator
  d=6 #conversion efficiency of biomass from prey to predator
  m=0.1 #predator mortality rate
  
  x0=0.675 #population density of prey at T=0
  y0=0.928 #population density of predator at T=0
  x=zeros(N,ss)
  y=zeros(N,ss)
  
  x[1,:].=x0
  y[1,:].=y0
  
  
  for j in 1:ss
    gen_rept[j]=sqrt(dt).*gen_rept[j]  
  end
  
  
  for j in 1:ss
    for i in 2:N
      x[i,j]=x[i-1,j].+(r.*x[i-1,j].*(1-x[i-1,j]).-((a.*x[i-1,j].*y[i-1,j])./(1 .+b.*x[i-1,j]))).*dt.+x[i-1,j].*sig_X.*gen_rept[j][i-1,1];
      y[i,j]=y[i-1,j].+(((d.*x[i-1,j].*y[i-1,j])./(1 .+b.*x[i-1,j])).*(1-y[i-1,j]./k1)-m.*y[i-1,j]).*dt.+y[i-1,j].*sig_Y.*gen_rept[j][i-1,2];
    end
  end
  
  return (std(x,dims=1),std(y,dims=1))
end
################# generate fluctuation size for constant k and varying q  ##############################

q=[0:0.1:1;] #values of species response correlation in sequence
k=0 # temporal environmental correlation
ss=1# Number of replicated time serise

deviate1=zeros(length(q),ss) #zero vector generated to store the  prey flactuation size for a given sequence of q
deviate2=zeros(length(q),ss) #zero vector generated to store the  predator flactuation size for a given sequence of q

for i in 1:length(q)
  standard_deviation=fluc(k,q[i],ss)
  deviate1[i,:]=standard_deviation[1]
  deviate2[i,:]=standard_deviation[2]
  println(i)
end


writedlm("mod1_sd_trend_q.txt",[deviate1 deviate2]) # save the data in text format
pp=readdlm("mod1_sd_trend_q.txt") #read the stored data


############ generate fluctuation size for constant q and varying k ######################################
q=0
k=[0:0.1:0.9;]
ss=1# Number of replicated time serise

deviate1=zeros(length(k),ss) #zero vector generated to store the  prey flactuation size for a given sequence of q
deviate2=zeros(length(k),ss) #zero vector generated to store the  predator flactuation size for a given sequence of q

for i in 1:length(k)
  standard_deviation=fluc(k[i],q,ss)
  deviate1[i,:]=standard_deviation[1]
  deviate2[i,:]=standard_deviation[2]
  println(i)
end


writedlm("mod1_sd_trend_k.txt",[deviate1 deviate2])
pp=readdlm("mod1_sd_trend_k.txt")


###########################################################################################
###########################################################################################
############  GENERATION OF FLUCTUATION SIZE FOR SHALLOW LAKE MODEL #######################
###########################################################################################
###########################################################################################

function fluc(k,q,ss)
  T=1500000
  dt=0.1
  N=convert(Int,T/dt)
  
  
  sig_X=0.02
  sig_Y=0.02
  
  
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
  
  
  r1=0.05 #maximum growth rate of Macrophyte
  r2=0.1  #maximum growth rate of Algae
  K=1     #carrying capacity ofthemacrophyte population
  h1=0.2  #half saturation constant on effect of algae on macrophyte
  h2=2    # half saturation constant on effect of macrophyte on algae
  p=4     # hill exponent
  T0=5.3  #Nutrient loading
  
  
  x0=0.7
  y0=4.408
  x=zeros(N,ss)
  y=zeros(N,ss)
  x[1,:].=x0
  y[1,:].=y0
  
  for j in 1:ss
    gen_rept[j]=sqrt(dt).*gen_rept[j]
  end
  
  
  for j in 1:ss
    for i in 2:N
      x[i,j]=x[i-1,j]+r1*x[i-1,j]*(1-(x[i-1,j]/(K*h2^p))*(h2^p+y[i-1,j]^p))*dt+x[i-1,j]*sig_X*gen_rept[j][i-1,1]
      y[i,j]=y[i-1,j]+r2*y[i-1,j]*(1-(y[i-1,j]/(h1*T0))*(h1+x[i-1,j]))*dt+y[i-1,j]*sig_Y*gen_rept[j][i-1,2]
    end
  end
  
  
  return (std(x,dims=1),std(y,dims=1))
end

########## generate fluctuation size for constant k and varying q  ##########

q=[0:0.1:1;] #values of species response correlation in sequence
k=0 # temporal environmental correlation
ss=1# Number of replicated time serise

deviate1=zeros(length(q),ss) #zero vector generated to store the  prey flactuation size for a given sequence of q
deviate2=zeros(length(q),ss) #zero vector generated to store the  predator flactuation size for a given sequence of q

for i in 1:length(q)
  standard_deviation=fluc(k,q[i],ss)
  deviate1[i,:]=standard_deviation[1]
  deviate2[i,:]=standard_deviation[2]
  println(i)
end


writedlm("mod2_sd_trend_q.txt",[deviate1 deviate2]) # save the data in text format
pp=readdlm("mod2_sd_trend_q.txt") #read the stored data


################# generate fluctuation size for constant q and varying k #############################
q=0
k=[0:0.1:0.9;]
ss=1# Number of replicated time serise

deviate1=zeros(length(k),ss) #zero vector generated to store the  prey flactuation size for a given sequence of q
deviate2=zeros(length(k),ss) #zero vector generated to store the  predator flactuation size for a given sequence of q

for i in 1:length(k)
  standard_deviation=fluc(k[i],q,ss)
  deviate1[i,:]=standard_deviation[1]
  deviate2[i,:]=standard_deviation[2]
  println(i)
end


writedlm("mod2_sd_trend_k.txt",[deviate1 deviate2])
pp=readdlm("mod2_sd_trend_k.txt")

###############################################################################################
###############################################################################################
############  GENERATION OF FLUCTUATION SIZE FOR MACROALGAE-CORAL MODEL #######################
###############################################################################################
###############################################################################################

function fluc(k,q,ss)
  T=3000000
  dt=0.1
  N=convert(Int,T/dt)
  
  
  sig_X=0.02
  sig_Y=0.02
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
  
  
  beta=0.35 #maximum grazing rate
  x0=0.0055
  y0=0.557
  
  x=zeros(N,ss)
  y=zeros(N,ss)
  x[1,:].=x0
  y[1,:].=y0

  
  for j in 1:ss
    gen_rept[j]=sqrt(dt).*gen_rept[j]
  end


  for j in 1:ss
    for i in 2:N
      x[i,j]=x[i-1,j].+(x[i-1,j].*(0.1.*y[i-1,j].-(beta./(1-y[i-1,j])).+0.77.*(1.0.-x[i-1,j].-y[i-1,j])).+0.005.*(1-x[i-1,j]-y[i-1,j])).*dt+x[i-1,j].*sig_X.*gen_rept[j][i-1,1]
      y[i,j]=y[i-1].+(y[i-1].*(0.55.*(1.0.-x[i-1].-y[i-1]).-0.24.-0.1.*x[i-1])).*dt.+y[i-1].*sig_Y.*gen_rept[j][i-1,2]
    end
  end

  return (std(x,dims=1),std(y,dims=1))
end

################# generate fluctuation size for constant k and varying q  ##############################
q=[0:0.1:1;]
k=0
ss=1# Number of replicated time serise

deviate1=zeros(length(q),ss) #zero vector generated to store the  prey flactuation size for a given sequence of q
deviate2=zeros(length(q),ss) #zero vector generated to store the  predator flactuation size for a given sequence of q

for i in 1:length(q)
  standard_deviation=fluc(k,q[i],ss)
  deviate1[i,:]=standard_deviation[1]
  deviate2[i,:]=standard_deviation[2]
  println(i)
end


writedlm("coral_sd_trend_q.txt",[deviate1 deviate2])
pp=readdlm("coral_sd_trend_q.txt")

############ generate fluctuation size for constant q and varying k ######################################

q=0
k=[0:0.1:0.9;]
ss=1# Number of replicated time serise

deviate1=zeros(length(k),ss) #zero vector generated to store the  prey flactuation size for a given sequence of q
deviate2=zeros(length(k),ss) #zero vector generated to store the  predator flactuation size for a given sequence of q

for i in 1:length(k)
  standard_deviation=fluc(k[i],q,ss)
  deviate1[i,:]=standard_deviation[1]
  deviate2[i,:]=standard_deviation[2]
  println(i)
end


writedlm("coral_sd_trend_k.txt",[deviate1 deviate2])
pp=readdlm("coral_sd_trend_k.txt")

################################################################################################################################
################################################################################################################################
############  GENERATION OF FLUCTUATION SIZE FOR TWO SPECIES MODEL WITH HOLLING TYPE FUNCTIONAL RESPONSE #######################
################################################################################################################################
################################################################################################################################

function fluc(k,q,ss)
  T=1500000
  dt=0.1
  N=convert(Int,T/dt)
  
  
  
  sig_X=0.03
  sig_Y=0.03
  
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
  r2=-0.5 
  a1=2
  a2=2   
  b1=2
  b2=2
  c1=1
  c2=1
  h1=0.5
  h2=0.5

  x0=1.66
  y0=1.99
  
  x=zeros(N,ss)
  y=zeros(N,ss)
  x[1,:].=x0
  y[1,:].=y0
  
  
  for j in 1:ss
    gen_rept[j]=sqrt(dt).*gen_rept[j]
  end
  
  for j in 1:ss
    for i in 2:N
      x[i,j]=x[i-1,j].+(x[i-1,j].*(r1.-c1.*x[i-1,j].+b1.*(a1.*y[i-1,j]./(1.0.+a1.*h1.*y[i-1,j])))).*dt.+x[i-1,j].*sig_X.*gen_rept[j][i-1,1]
      y[i,j]=y[i-1,j].+(y[i-1,j].*(r2.-c2.*y[i-1,j].+b2.*(a2.*x[i-1,j]./(1.0.+a2.*h2.*x[i-1,j])))).*dt.+y[i-1,j].*sig_Y.*gen_rept[j][i-1,2]
    end
  end

  return (std(x,dims=1),std(y,dims=1))
end

################# generate fluctuation size for constant k and varying q  ##############################
q=[0:0.1:1;]
k=0
ss=1# Number of replicated time serise

deviate1=zeros(length(q),ss) #zero vector generated to store the  prey flactuation size for a given sequence of q
deviate2=zeros(length(q),ss) #zero vector generated to store the  predator flactuation size for a given sequence of q

for i in 1:length(q)
  standard_deviation=fluc(k,q[i],ss)
  deviate1[i,:]=standard_deviation[1]
  deviate2[i,:]=standard_deviation[2]
  println(i)
end


writedlm("mod4_ob_sd_trend_q.txt",[deviate1 deviate2])
pp=readdlm("mod4_ob_sd_trend_q.txt")

############ generate fluctuation size for constant q and varying k ######################################

q=0
k=[0:0.1:0.9;]
ss=1# Number of replicated time serise

deviate1=zeros(length(k),ss) #zero vector generated to store the  prey flactuation size for a given sequence of q
deviate2=zeros(length(k),ss) #zero vector generated to store the  predator flactuation size for a given sequence of q

for i in 1:length(k)
  standard_deviation=fluc(k[i],q,ss)
  deviate1[i,:]=standard_deviation[1]
  deviate2[i,:]=standard_deviation[2]
  println(i)
end


writedlm("mod4_ob_sd_trend_k.txt",[deviate1 deviate2])
pp=readdlm("mod4_ob_sd_trend_k.txt")

###########################################################################################################
###########################################################################################################
############  GENERATION OF FLUCTUATION SIZE FOR PLANT-POLLINATOR INTERACTION MODEL #######################
###########################################################################################################
###########################################################################################################
function fluc(k,q,ss)
  T=2500000
  dt=0.1
  N=convert(Int,T/dt)
  
  
  
  sig_X=0.005
  sig_Y=0.005
  
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
  
  
  a0=2.06 #per-plant attack rate on rewards
  x0=0.1519
  y0=0.2562
  
  x=zeros(N,ss)
  y=zeros(N,ss)
  x[1,:].=x0
  y[1,:].=y0
  
  
  for j in 1:ss
    gen_rept[j]=sqrt(dt).*gen_rept[j]
  end
  
  for j in 1:ss
    for i in 2:N
      x[i,j]=x[i-1,j].+(x[i-1,j].*(1.42.*(0.5.+0.5.*(a0.*x[i-1,j].*y[i-1,j]./(1.0.+a0.*0.01.*x[i-1,j].+a0.*y[i-1,j].*x[i-1,j]))).-0.61.*x[i-1,j].-0.67)).*dt.+x[i-1,j].*sig_X.*gen_rept[j][i-1,1]
      y[i,j]=y[i-1,j].+(y[i-1,j].*(1.0.+a0.*x[i-1,j]./(1.0.+a0.*0.01.*x[i-1,j]).-2.0.*y[i-1,j].-0.8)).*dt.+y[i-1,j].*sig_Y.*gen_rept[j][i-1,2]
    end
  end

  return (std(x,dims=1),std(y,dims=1))
end

################# generate fluctuation size for constant k and varying q  ##############################

q=[0:0.1:1;]
k=0
ss=1# Number of replicated time serise

deviate1=zeros(length(q),ss) #zero vector generated to store the  prey flactuation size for a given sequence of q
deviate2=zeros(length(q),ss) #zero vector generated to store the  predator flactuation size for a given sequence of q

for i in 1:length(q)
  standard_deviation=fluc(k,q[i],ss)
  deviate1[i,:]=standard_deviation[1]
  deviate2[i,:]=standard_deviation[2]
  println(i)
end

writedlm("mod5_bf_sd_trend_q.txt",[deviate1 deviate2])
pp=readdlm("mod5_bf_sd_trend_q.txt")

############ generate fluctuation size for constant q and varying k ######################################

q=0
k=[0:0.1:0.9;]
ss=1# Number of replicated time serise

deviate1=zeros(length(k),ss) #zero vector generated to store the  prey flactuation size for a given sequence of q
deviate2=zeros(length(k),ss) #zero vector generated to store the  predator flactuation size for a given sequence of q

for i in 1:length(k)
  standard_deviation=fluc(k[i],q,ss)
  deviate1[i,:]=standard_deviation[1]
  deviate2[i,:]=standard_deviation[2]
  println(i)
end


writedlm("mod5_bf_sd_trend_k.txt",[deviate1 deviate2])
pp=readdlm("mod5_bf_sd_trend_k.txt")

###########################################################################################################
###########################################################################################################
############  GENERATION OF FLUCTUATION SIZE FOR PATCH OCCUPANCY MODEL OF MUTUALISM #######################
###########################################################################################################
###########################################################################################################

function fluc(k,q,ss)
  T=1500000
  dt=0.1
  N=convert(Int,T/dt)

  sig_X=0.03
  sig_Y=0.03
  
  
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
  
  
  d0=0.2 #fraction of destroyed habitat
  x0=0.72
  y0=0.47
  
  x=zeros(N,ss)
  y=zeros(N,ss)

  x[1,:].=x0
  y[1,:].=y0
  
  
  for j in 1:ss
    gen_rept[j]=sqrt(dt).*gen_rept[j]
  end
  
  for j in 1:ss
    for i in 2:N
      x[i,j]=x[i-1,j]+(2.0*y[i-1,j]*(1-d0-x[i-1,j])-0.1*x[i-1,j])*dt+x[i-1,j]*sig_X*gen_rept[j][i-1,1]
      y[i,j]=y[i-1,j]+(2.0*y[i-1,j]*(x[i-1,j]-y[i-1,j])-0.5*y[i-1,j])*dt+y[i-1,j]*sig_Y*gen_rept[j][i-1,2]
    end
  end
  
  return (std(x,dims=1),std(y,dims=1))
end


################# generate fluctuation size for constant k and varying q  ##############################

q=[0:0.1:1;]
k=0
ss=1# Number of replicated time serise

deviate1=zeros(length(q),ss) #zero vector generated to store the  prey flactuation size for a given sequence of q
deviate2=zeros(length(q),ss) #zero vector generated to store the  predator flactuation size for a given sequence of q

for i in 1:length(q)
  standard_deviation=fluc(k,q[i],ss)
  deviate1[i,:]=standard_deviation[1]
  deviate2[i,:]=standard_deviation[2]
  println(i)
end


writedlm("mod6_spmut_sd_trend_q.txt",[deviate1 deviate2])
pp=readdlm("mod6_spmut_sd_trend_q.txt")

############ generate fluctuation size for constant q and varying k ######################################

q=0
k=[0:0.1:0.9;]
ss=1# Number of replicated time serise

deviate1=zeros(length(k),ss) #zero vector generated to store the  prey flactuation size for a given sequence of q
deviate2=zeros(length(k),ss) #zero vector generated to store the  predator flactuation size for a given sequence of q

for i in 1:length(k)
  standard_deviation=fluc(k[i],q,ss)
  deviate1[i,:]=standard_deviation[1]
  deviate2[i,:]=standard_deviation[2]
  println(i)
end


writedlm("mod6_spmut_sd_trend_k.txt",[deviate1 deviate2])
pp=readdlm("mod6_spmut_sd_trend_k.txt")

