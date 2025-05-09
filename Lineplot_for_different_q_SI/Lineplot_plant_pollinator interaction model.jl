###################### Line plot Figure.C2(b) in SI #######################################
################### For plant pollinator interaction model#######################
######################################################################################


using Distributed
nprocs()
addprocs(15)
w=workers()

using Plots
using DataFrames
using Statistics
using StatsPlots
using VectorizedStatistics
using LaTeXStrings
using DelimitedFiles
using CSV

@everywhere using Statistics
@everywhere function gen(q,ss)

  T=300000 #500000
  dt=0.1
  N=convert(Int,T/dt)
  
  sig_X=0.006
  sig_Y=0.006
  k=0.8
  
  

  
  sub_mat=zeros(N,ss)

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


    a0=2.0
   
  
    
    x0=0.14
    y0=0.24


    x=zeros(N)
    y=zeros(N)
    a=zeros(N)



    x[1]=x0
    y[1]=y0
    a[1]=a0


    noise.=sqrt(dt).*noise
      
      
    
    for i in 2:N
        x[i]=x[i-1].+(x[i-1].*(1.42.*(0.5.+0.5.*(a[i-1].*x[i-1].*y[i-1]./(1.0.+a[i-1].*0.01.*x[i-1].+a[i-1].*y[i-1].*x[i-1]))).-0.61.*x[i-1].-0.67)).*dt.+x[i-1].*sig_X.*noise[i-1,1]
        y[i]=y[i-1].+(y[i-1].*(1.0.+a[i-1].*x[i-1]./(1.0.+a[i-1].*0.01.*x[i-1]).-2.0.*y[i-1].-0.8)).*dt.+y[i-1].*sig_Y.*noise[i-1,2]
        a[i]=a[i-1].+((2.12-2.0)/T)*dt
    end

    
    sub_mat[:,ll].=x
    println(ll)
  end

  return sub_mat
  
end


T=300000#500000
dt=0.1
N=convert(Int,T/dt)
m=zeros(N)
m[1]=2.0
for i in 2:N
  m[i]=m[i-1].+((2.12-2.0)/T)*dt
end

rep=100
q=0.0


t1=@spawnat w[1] gen(q,rep)
t2=@spawnat w[2] gen(q,rep)
t3=@spawnat w[3] gen(q,rep)
t4=@spawnat w[4] gen(q,rep)
t5=@spawnat w[5] gen(q,rep)
t6=@spawnat w[6] gen(q,rep)
t7=@spawnat w[7] gen(q,rep)
t8=@spawnat w[8] gen(q,rep)
t9=@spawnat w[9] gen(q,rep)
t10=@spawnat w[10] gen(q,rep)
t11=@spawnat w[11] gen(q,rep)
t12=@spawnat w[12] gen(q,rep)
t13=@spawnat w[13] gen(q,rep)
t14=@spawnat w[14] gen(q,rep)
t15=@spawnat w[15] gen(q,rep)

sub1=fetch(t1)
sub2=fetch(t2)
sub3=fetch(t3)
sub4=fetch(t4)
sub5=fetch(t5)
sub6=fetch(t6)
sub7=fetch(t7)
sub8=fetch(t8)
sub9=fetch(t9)
sub10=fetch(t10)
sub11=fetch(t11)
sub12=fetch(t12)
sub13=fetch(t13)
sub14=fetch(t14)
sub15=fetch(t15)

mat=hcat(sub1,sub2,sub3,sub4,sub5,sub6,sub7,sub8,sub9,sub10,sub11,sub12,sub13,sub14,sub15)



M=mat[1:100:end,:]
R=m[1:100:end,:]

mn=mean(M,dims=2)
p5 = vquantile(copy(M),0.05,dims=2)  
p95 = vquantile(copy(M), 0.95; dims=2)  

data=DataFrame(para=vec(R[:,1]),mn=vec(mn),p5s=vec(p5),p95s=vec(p95))


@df data plot(:para,:mn,ribbon=(:mn-:p5s,-:mn+:p95s),fillalpha=0.1,framestyle=:box,label="q=0.0",foreground_color_legend=nothing,
legendfontsize=11,linewidth=1.75,thickness_scaling=1.1,legend=:outertop,legend_column=-1,
grid=false,xlab=L"Per-plant\; attack\; rate\;on\;rewards\; (a)",widen=false,color=:red2,
tickfontsize=13,guidefontsize=17,ylab=L"Plant\;abundance\; (P)",xlim=(2,2.12),ylim=(0.12,0.75))




q=0.45

t1=@spawnat w[1] gen(q,rep)
t2=@spawnat w[2] gen(q,rep)
t3=@spawnat w[3] gen(q,rep)
t4=@spawnat w[4] gen(q,rep)
t5=@spawnat w[5] gen(q,rep)
t6=@spawnat w[6] gen(q,rep)
t7=@spawnat w[7] gen(q,rep)
t8=@spawnat w[8] gen(q,rep)
t9=@spawnat w[9] gen(q,rep)
t10=@spawnat w[10] gen(q,rep)
t11=@spawnat w[11] gen(q,rep)
t12=@spawnat w[12] gen(q,rep)
t13=@spawnat w[13] gen(q,rep)
t14=@spawnat w[14] gen(q,rep)
t15=@spawnat w[15] gen(q,rep)

sub1=fetch(t1)
sub2=fetch(t2)
sub3=fetch(t3)
sub4=fetch(t4)
sub5=fetch(t5)
sub6=fetch(t6)
sub7=fetch(t7)
sub8=fetch(t8)
sub9=fetch(t9)
sub10=fetch(t10)
sub11=fetch(t11)
sub12=fetch(t12)
sub13=fetch(t13)
sub14=fetch(t14)
sub15=fetch(t15)

mat=hcat(sub1,sub2,sub3,sub4,sub5,sub6,sub7,sub8,sub9,sub10,sub11,sub12,sub13,sub14,sub15)


M=mat[1:100:end,:]
R=m[1:100:end,:]

mn=mean(M,dims=2)
p5 = vquantile(copy(M),0.05,dims=2)  
p95 = vquantile(copy(M), 0.95; dims=2)  

data=DataFrame(para=vec(R[:,1]),mn=vec(mn),p5s=vec(p5),p95s=vec(p95))



@df data plot!(:para,:mn,ribbon=(:mn-:p5s,-:mn+:p95s),fillalpha=0.1,framestyle=:box,label="q=0.45",foreground_color_legend=nothing,
legendfontsize=11,linewidth=1.75,thickness_scaling=1.1,legend=:outertop,legend_column=-1,
grid=false,xlab=L"Per-plant\; attack\; rate\;on\;rewards\; (a)",widen=false,color=:darkorange,
tickfontsize=13,guidefontsize=17,ylab=L"Plant\;abundance\; (P)",xlim=(2,2.12),ylim=(0.12,0.75))

q=0.9

t1=@spawnat w[1] gen(q,rep)
t2=@spawnat w[2] gen(q,rep)
t3=@spawnat w[3] gen(q,rep)
t4=@spawnat w[4] gen(q,rep)
t5=@spawnat w[5] gen(q,rep)
t6=@spawnat w[6] gen(q,rep)
t7=@spawnat w[7] gen(q,rep)
t8=@spawnat w[8] gen(q,rep)
t9=@spawnat w[9] gen(q,rep)
t10=@spawnat w[10] gen(q,rep)
t11=@spawnat w[11] gen(q,rep)
t12=@spawnat w[12] gen(q,rep)
t13=@spawnat w[13] gen(q,rep)
t14=@spawnat w[14] gen(q,rep)
t15=@spawnat w[15] gen(q,rep)

sub1=fetch(t1)
sub2=fetch(t2)
sub3=fetch(t3)
sub4=fetch(t4)
sub5=fetch(t5)
sub6=fetch(t6)
sub7=fetch(t7)
sub8=fetch(t8)
sub9=fetch(t9)
sub10=fetch(t10)
sub11=fetch(t11)
sub12=fetch(t12)
sub13=fetch(t13)
sub14=fetch(t14)
sub15=fetch(t15)

mat=hcat(sub1,sub2,sub3,sub4,sub5,sub6,sub7,sub8,sub9,sub10,sub11,sub12,sub13,sub14,sub15)


M=mat[1:100:end,:]
R=m[1:100:end,:]

mn=mean(M,dims=2)
p5 = vquantile(copy(M),0.05,dims=2)  
p95 = vquantile(copy(M), 0.95; dims=2)  

data=DataFrame(para=vec(R[:,1]),mn=vec(mn),p5s=vec(p5),p95s=vec(p95))


@df data plot!(:para,:mn,ribbon=(:mn-:p5s,-:mn+:p95s),fillalpha=0.1,framestyle=:box,label="q=0.9",foreground_color_legend=nothing,
legendfontsize=11,linewidth=1.75,thickness_scaling=1.1,legend=:outertop,legend_column=-1,
grid=false,xlab=L"Per-plant\; attack\; rate\;on\;rewards\; (a)",widen=false,color=:darkgreen,
tickfontsize=13,guidefontsize=17,ylab=L"Plant\;abundance\; (P)",xlim=(2,2.12),ylim=(0.12,0.75))


vline!([2.112],linewidth=1.5,linestyle=:dash,color=:darkgrey,label="",dpi=750)


## For reverse direction of the parameter
@everywhere using Statistics
@everywhere function gen_rev(q,ss)

  T=300000#500000
  dt=0.1
  N=convert(Int,T/dt)
  
  sig_X=0.006
  sig_Y=0.006
  k=0.8
  
  

  
  sub_mat=zeros(N,ss)

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


    a0=2.12
   
    x0=0.716
    y0=0.848
    


    x=zeros(N)
    y=zeros(N)
    a=zeros(N)



    x[1]=x0
    y[1]=y0
    a[1]=a0


    noise.=sqrt(dt).*noise
      
      
    
    for i in 2:N
        x[i]=x[i-1].+(x[i-1].*(1.42.*(0.5.+0.5.*(a[i-1].*x[i-1].*y[i-1]./(1.0.+a[i-1].*0.01.*x[i-1].+a[i-1].*y[i-1].*x[i-1]))).-0.61.*x[i-1].-0.67)).*dt.+x[i-1].*sig_X.*noise[i-1,1]
        y[i]=y[i-1].+(y[i-1].*(1.0.+a[i-1].*x[i-1]./(1.0.+a[i-1].*0.01.*x[i-1]).-2.0.*y[i-1].-0.8)).*dt.+y[i-1].*sig_Y.*noise[i-1,2]
        a[i]=a[i-1].-((2.12-2.0)/T)*dt
    end

    
    sub_mat[:,ll].=x
    println(ll)
  end

  return sub_mat
  
end


T=300000#500000
dt=0.1
N=convert(Int,T/dt)
m=zeros(N)
m[1]=2.12
for i in 2:N
  m[i]=m[i-1].-((2.12-2.0)/T)*dt
end

rep=100
q=0.0



t1=@spawnat w[1] gen_rev(q,rep)
t2=@spawnat w[2] gen_rev(q,rep)
t3=@spawnat w[3] gen_rev(q,rep)
t4=@spawnat w[4] gen_rev(q,rep)
t5=@spawnat w[5] gen_rev(q,rep)
t6=@spawnat w[6] gen_rev(q,rep)
t7=@spawnat w[7] gen_rev(q,rep)
t8=@spawnat w[8] gen_rev(q,rep)
t9=@spawnat w[9] gen_rev(q,rep)
t10=@spawnat w[10] gen_rev(q,rep)
t11=@spawnat w[11] gen_rev(q,rep)
t12=@spawnat w[12] gen_rev(q,rep)
t13=@spawnat w[13] gen_rev(q,rep)
t14=@spawnat w[14] gen_rev(q,rep)
t15=@spawnat w[15] gen_rev(q,rep)

sub1=fetch(t1)
sub2=fetch(t2)
sub3=fetch(t3)
sub4=fetch(t4)
sub5=fetch(t5)
sub6=fetch(t6)
sub7=fetch(t7)
sub8=fetch(t8)
sub9=fetch(t9)
sub10=fetch(t10)
sub11=fetch(t11)
sub12=fetch(t12)
sub13=fetch(t13)
sub14=fetch(t14)
sub15=fetch(t15)

mat=hcat(sub1,sub2,sub3,sub4,sub5,sub6,sub7,sub8,sub9,sub10,sub11,sub12,sub13,sub14,sub15)



M=mat[1:100:end,:]
R=m[1:100:end,:]

mn=mean(M,dims=2)
p5 = vquantile(copy(M),0.05,dims=2)  
p95 = vquantile(copy(M), 0.95; dims=2)  

data=DataFrame(para=vec(R[:,1]),mn=vec(mn),p5s=vec(p5),p95s=vec(p95))


@df data plot!(:para,:mn,ribbon=(:mn-:p5s,-:mn+:p95s),fillalpha=0.1,framestyle=:box,label="",foreground_color_legend=nothing,
legendfontsize=11,linewidth=1.75,thickness_scaling=1.1,legend=:outertop,legend_column=-1,
grid=false,xlab=L"Per-plant\; attack\; rate\;on\;rewards\; (a)",widen=false,color=:red2,
tickfontsize=13,guidefontsize=17,ylab=L"Plant\;abundance\; (P)",xlim=(2,2.12),ylim=(0.12,0.75))




q=0.45

t1=@spawnat w[1] gen_rev(q,rep)
t2=@spawnat w[2] gen_rev(q,rep)
t3=@spawnat w[3] gen_rev(q,rep)
t4=@spawnat w[4] gen_rev(q,rep)
t5=@spawnat w[5] gen_rev(q,rep)
t6=@spawnat w[6] gen_rev(q,rep)
t7=@spawnat w[7] gen_rev(q,rep)
t8=@spawnat w[8] gen_rev(q,rep)
t9=@spawnat w[9] gen_rev(q,rep)
t10=@spawnat w[10] gen_rev(q,rep)
t11=@spawnat w[11] gen_rev(q,rep)
t12=@spawnat w[12] gen_rev(q,rep)
t13=@spawnat w[13] gen_rev(q,rep)
t14=@spawnat w[14] gen_rev(q,rep)
t15=@spawnat w[15] gen_rev(q,rep)

sub1=fetch(t1)
sub2=fetch(t2)
sub3=fetch(t3)
sub4=fetch(t4)
sub5=fetch(t5)
sub6=fetch(t6)
sub7=fetch(t7)
sub8=fetch(t8)
sub9=fetch(t9)
sub10=fetch(t10)
sub11=fetch(t11)
sub12=fetch(t12)
sub13=fetch(t13)
sub14=fetch(t14)
sub15=fetch(t15)

mat=hcat(sub1,sub2,sub3,sub4,sub5,sub6,sub7,sub8,sub9,sub10,sub11,sub12,sub13,sub14,sub15)


M=mat[1:100:end,:]
R=m[1:100:end,:]

mn=mean(M,dims=2)
p5 = vquantile(copy(M),0.05,dims=2)  
p95 = vquantile(copy(M), 0.95; dims=2)  

data=DataFrame(para=vec(R[:,1]),mn=vec(mn),p5s=vec(p5),p95s=vec(p95))


@df data plot!(:para,:mn,ribbon=(:mn-:p5s,-:mn+:p95s),fillalpha=0.1,framestyle=:box,label="",foreground_color_legend=nothing,
legendfontsize=11,linewidth=1.75,thickness_scaling=1.1,legend=:outertop,legend_column=-1,
grid=false,xlab=L"Per-plant\; attack\; rate\;on\;rewards\; (a)",widen=false,color=:darkorange,
tickfontsize=13,guidefontsize=17,ylab=L"Plant\;abundance\; (P)",xlim=(2,2.12),ylim=(0.12,0.75))

q=0.9

t1=@spawnat w[1] gen_rev(q,rep)
t2=@spawnat w[2] gen_rev(q,rep)
t3=@spawnat w[3] gen_rev(q,rep)
t4=@spawnat w[4] gen_rev(q,rep)
t5=@spawnat w[5] gen_rev(q,rep)
t6=@spawnat w[6] gen_rev(q,rep)
t7=@spawnat w[7] gen_rev(q,rep)
t8=@spawnat w[8] gen_rev(q,rep)
t9=@spawnat w[9] gen_rev(q,rep)
t10=@spawnat w[10] gen_rev(q,rep)
t11=@spawnat w[11] gen_rev(q,rep)
t12=@spawnat w[12] gen_rev(q,rep)
t13=@spawnat w[13] gen_rev(q,rep)
t14=@spawnat w[14] gen_rev(q,rep)
t15=@spawnat w[15] gen_rev(q,rep)

sub1=fetch(t1)
sub2=fetch(t2)
sub3=fetch(t3)
sub4=fetch(t4)
sub5=fetch(t5)
sub6=fetch(t6)
sub7=fetch(t7)
sub8=fetch(t8)
sub9=fetch(t9)
sub10=fetch(t10)
sub11=fetch(t11)
sub12=fetch(t12)
sub13=fetch(t13)
sub14=fetch(t14)
sub15=fetch(t15)

mat=hcat(sub1,sub2,sub3,sub4,sub5,sub6,sub7,sub8,sub9,sub10,sub11,sub12,sub13,sub14,sub15)


M=mat[1:100:end,:]
R=m[1:100:end,:]

mn=mean(M,dims=2)
p5 = vquantile(copy(M),0.05,dims=2)  
p95 = vquantile(copy(M), 0.95; dims=2)  

data=DataFrame(para=vec(R[:,1]),mn=vec(mn),p5s=vec(p5),p95s=vec(p95))


@df data plot!(:para,:mn,ribbon=(:mn-:p5s,-:mn+:p95s),fillalpha=0.1,framestyle=:box,label="",foreground_color_legend=nothing,
legendfontsize=11,linewidth=1.75,thickness_scaling=1.1,legend=:outertop,legend_column=-1,
grid=false,xlab=L"Per-plant\; attack\; rate\;on\;rewards\; (a)",widen=false,color=:darkgreen,
tickfontsize=13,guidefontsize=17,ylab=L"Plant\;abundance\; (P)",xlim=(2,2.12),ylim=(0.12,0.75))


vline!([2.011],linewidth=1.5,linestyle=:dash,color=:darkgrey,label="",dpi=750)

#savefig("mut_oldmod_95shaded_comb.png")




