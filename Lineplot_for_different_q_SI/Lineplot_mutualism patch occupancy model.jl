###################### Line plot Figure.C2(e) in SI #######################################
################### For plant-animal patch occupancy mutualism model#######################
###########################################################################################


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
    
    T=200000
    dt=0.1
    N=convert(Int,T/dt)
    sig_X=0.03
    sig_Y=0.03
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


        m0=0.0
        x0=0.931
        y0=0.681



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
            m[i]=m[i-1]+(0.5/T)*dt
        end
        
        sub_mat[:,ll].=x
        println(ll)
    end
    
    return sub_mat
  
end

T=200000
dt=0.1
N=convert(Int,T/dt)
m=zeros(N)
m[1]=0.0
for i in 2:N
  m[i]=m[i-1].+((0.5)/T)*dt
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



M=mat[1:500:end,:]
R=m[1:500:end,:]

mn=mean(M,dims=2)
p5 = vquantile(copy(M),0.05,dims=2)  
p95 = vquantile(copy(M), 0.95; dims=2)  

data=DataFrame(para=vec(R[:,1]),mn=vec(mn),p5s=vec(p5),p95s=vec(p95))



@df data plot(:para,:mn,ribbon=(:mn-:p5s,-:mn+:p95s),fillalpha=0.1,framestyle=:box,label="q=0.0",foreground_color_legend=nothing,
legendfontsize=11,linewidth=2.5,thickness_scaling=1.1,legend=:outertop,legend_column=-1,
grid=false,xlab=L"Fraction\; of \; destroyrd\;habitat\; (d)",widen=false,color=:red2,
tickfontsize=13,guidefontsize=17,ylab=L"Fraction\; of \;habitat\; occupied\; by\; Plant\;(p)",xlim=(0.2,0.5),ylim=(0.0,0.8))




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



M=mat[1:500:end,:]
R=m[1:500:end,:]

mn=mean(M,dims=2)
p5 = vquantile(copy(M),0.05,dims=2)  
p95 = vquantile(copy(M), 0.95; dims=2)  

data=DataFrame(para=vec(R[:,1]),mn=vec(mn),p5s=vec(p5),p95s=vec(p95))


@df data plot!(:para,:mn,ribbon=(:mn-:p5s,-:mn+:p95s),fillalpha=0.1,framestyle=:box,label="q=0.45",foreground_color_legend=nothing,
legendfontsize=11,linewidth=2.5,thickness_scaling=1.1,legend=:outertop,legend_column=-1,
grid=false,xlab=L"Fraction\; of \; destroyrd\;habitat\; (d)",widen=false,color=:darkorange,
tickfontsize=13,guidefontsize=17,ylab=L"Fraction\; of \;habitat\; occupied\; by\; Plant\;(p)",xlim=(0.2,0.5),ylim=(0.0,0.8))



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



M=mat[1:500:end,:]
R=m[1:500:end,:]

mn=mean(M,dims=2)
p5 = vquantile(copy(M),0.05,dims=2)  
p95 = vquantile(copy(M), 0.95; dims=2)  

data=DataFrame(para=vec(R[:,1]),mn=vec(mn),p5s=vec(p5),p95s=vec(p95))

@df data plot!(:para,:mn,ribbon=(:mn-:p5s,-:mn+:p95s),fillalpha=0.1,framestyle=:box,label="q=0.9",foreground_color_legend=nothing,
legendfontsize=11,linewidth=2.5,thickness_scaling=1.1,legend=:outertop,legend_column=-1,
grid=false,xlab=L"Fraction\; of \; destroyed\;habitat\; (d)",widen=false,color=:darkgreen,
tickfontsize=13,guidefontsize=13,ylab=L"Fraction\; of \;habitat\; occupied\; by\; Plants\;(p)",xlim=(0.2,0.5),ylim=(0.0,0.8))



vline!([0.476],linewidth=1.5,linestyle=:dash,color=:darkgrey,label="",dpi=750)


#savefig("mut_spatial_mod6_combined.png")