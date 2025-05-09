
using LaTeXStrings

using Plots, Plots.PlotMeasures; gr()
using DelimitedFiles

function twiny(sp::Plots.Subplot)
    sp[:top_margin] = max(sp[:top_margin], 30Plots.px)
    plot!(sp.plt, inset = (sp[:subplot_index], bbox(0,0,1,1)))
    twinsp = sp.plt.subplots[end]
    twinsp[:xaxis][:mirror] = true
    twinsp[:background_color_inside] = RGBA{Float64}(0,0,0,0)
    Plots.link_axes!(sp[:yaxis], twinsp[:yaxis])
    twinsp
end
twiny(plt::Plots.Plot = current()) = twiny(plt[1])




###############################----------Plot flactuation size for Prey-Predator model----------------######################### 


####################################------Plot flactuation size w.r.t species response correlation------###################
pp=readdlm("mod1_sd_trend_q.txt")
q=[0:0.1:1;]
k=0

plot([pp[1,1]],label=L"\textbf{y}\;\;\;",c=:red,lw=5,foreground_color_legend=nothing,legendfontsize=15,legend=:outertop)
plot!(q,pp[:,1],linewidth=5,box=:on,widen=false,grid=false,label=L"\textbf{x}\;",legendfontsize=15,legend=:outertop,legend_column=-1,
tickfontsize=12,xlabel=L"Species\; response \;correlation \;(q)",foreground_color_legend=nothing,
ylab=L"Fluctuation\; size",guidefontsize=18,yaxis=(formatter=y->string(round( y / 10^-2;digits=3))),c=:blue,
y_foreground_color_axis=:blue,y_foreground_color_border=:blue,xlim=(0,1.0))
annotate!([(0.0,maximum(pp[:,1])*1.008,text(L"\times10^{-2}", 14, :black))])
annotate!([(1.04,maximum(pp[:,1])*1.008,text(L"\times10^{-2}", 14, :black))])

plot!(twinx(),q,pp[:,2],linewidth=5,color=:red,widen=false,grid=false,xlim=(0,1.0),
yaxis=(formatter=y->string(round( y / 10^-2;digits=3))),legendfontsize=11,legend=:none,
tickfontsize=12,y_foreground_color_axis=:red,y_foreground_color_border=:red,label="",dpi=750)

plot!(twiny(),xticks=:false,yticks=:false,y_foreground_color_border=:blue)
plot!(legend=(0.45,1.05),top_margin=10mm)

savefig("trend_std_m1_q.png")


####################################------Plot flactuation size w.r.t tempral environmental correlation------###################
pp=readdlm("mod1_sd_trend_k.txt")
q=0
k=[0:0.1:0.9;]

plot([pp[1,1]],label=L"\textbf{y}\;\;\;",c=:red,lw=5,foreground_color_legend=nothing,legendfontsize=15,legend=:outertop)
plot!(k,pp[:,1],linewidth=5,box=:on,widen=false,grid=false,label=L"\textbf{x}\;",legendfontsize=15,legend=:outertop,legend_column=-1,
tickfontsize=12,xlabel=L"Temporal \;correlation \;(k)",foreground_color_legend=nothing,
ylab=L"Fluctuation\; size",guidefontsize=18,yaxis=(formatter=y->string(round( y / 10^-1;digits=3))),c=:blue,
y_foreground_color_axis=:blue,y_foreground_color_border=:blue,xlim=(0,0.9))
annotate!([(0.0,maximum(pp[:,1])*1.05,text(L"\times10^{-1}", 14, :black))])
annotate!([(0.93,maximum(pp[:,1])*1.05,text(L"\times10^{-1}", 14, :black))])

plot!(twinx(),k,pp[:,2],linewidth=5,color=:red,widen=false,grid=false,xlim=(0,0.9),
yaxis=(formatter=y->string(round( y / 10^-1;digits=3))),legendfontsize=11,legend=:none,
tickfontsize=12,y_foreground_color_axis=:red,y_foreground_color_border=:red,label="",dpi=750)

plot!(twiny(),xticks=:false,yticks=:false,y_foreground_color_border=:blue)
plot!(legend=(0.45,1.05),top_margin=10mm)

savefig("trend_std_m1_k.png")

###############################----------Plot flactuation size for Shallow lake model----------------######################### 


####################################------Plot flactuation size w.r.t species response correlation------###################
q=[0:0.1:1;]
k=0
pp=readdlm("mod2_sd_trend_q.txt")

plot([pp[1,1]],label=L"\textbf{A}\;\;\;",c=:red,lw=5,foreground_color_legend=nothing,legendfontsize=15,legend=:outertop)
plot!(q,pp[:,1],linewidth=5,box=:on,widen=false,grid=false,label=L"\textbf{M}\;",legendfontsize=15,legend=:outertop,legend_column=-1,
tickfontsize=12,xlabel=L"Species\; response \;correlation \;(q)",foreground_color_legend=nothing,
ylab=L"Fluctuation\; size",guidefontsize=18,yaxis=(formatter=y->string(round( y / 10^-2;digits=3))),c=:blue,
y_foreground_color_axis=:blue,y_foreground_color_border=:blue,xlim=(0,1.0))
annotate!([(0.0,maximum(pp[:,1])*1.004,text(L"\times10^{-2}", 14, :black))])
annotate!([(1.04,maximum(pp[:,1])*1.004,text(L"\times10^{-2}", 14, :black))])

plot!(twinx(),q,pp[:,2],linewidth=5,color=:red,widen=false,grid=false,xlim=(0,1.0),
yaxis=(formatter=y->string(round( y / 10^-2;digits=3))),legendfontsize=11,legend=:none,
tickfontsize=12,y_foreground_color_axis=:red,y_foreground_color_border=:red,label="",dpi=750)

plot!(twiny(),xticks=:false,yticks=:false,y_foreground_color_border=:blue)
plot!(legend=(0.45,1.05),top_margin=10mm)

savefig("trend_std_m2_q.png")


####################################------Plot flactuation size w.r.t tempral environmental correlation------###################

pp=readdlm("mod2_sd_trend_k.txt")
q=0
k=[0:0.1:0.9;]

plot([pp[1,1]],label=L"\textbf{A}\;\;\;",c=:red,lw=5,foreground_color_legend=nothing,legendfontsize=15,legend=:outertop)
plot!(k,pp[:,1],linewidth=5,box=:on,widen=false,grid=false,label=L"\textbf{M}\;",legendfontsize=15,legend=:outertop,legend_column=-1,
tickfontsize=12,xlabel=L"Temporal \;correlation \;(k)",foreground_color_legend=nothing,
ylab=L"Fluctuation\; size",guidefontsize=18,yaxis=(formatter=y->string(round( y / 10^-1;digits=3))),c=:blue,
y_foreground_color_axis=:blue,y_foreground_color_border=:blue,xlim=(0,0.9))
annotate!([(0.0,maximum(pp[:,1])*1.05,text(L"\times10^{-1}", 14, :black))])
annotate!([(0.93,maximum(pp[:,1])*1.05,text(L"\times10^{-1}", 14, :black))])

plot!(twinx(),k,pp[:,2],linewidth=5,color=:red,widen=false,grid=false,xlim=(0,0.9),
yaxis=(formatter=y->string(round( y / 10^-1;digits=3))),legendfontsize=11,legend=:none,
tickfontsize=12,y_foreground_color_axis=:red,y_foreground_color_border=:red,label="",dpi=750)

plot!(twiny(),xticks=:false,yticks=:false,y_foreground_color_border=:blue)
plot!(legend=(0.45,1.05),top_margin=10mm)

savefig("trend_std_m2_k.png")


###############################----------Plot flactuation size for Coral-Macroalgae model----------------######################### 


####################################------Plot flactuation size w.r.t species response correlation------###################

q=[0:0.1:1;]
k=0
pp=readdlm("coral_sd_trend_q.txt")

plot([pp[1,1]],label=L"\textbf{C}\;\;\;",c=:red,lw=5,foreground_color_legend=nothing,legendfontsize=15,legend=:outertop)
plot!(q,pp[:,1],linewidth=5,box=:on,widen=false,grid=false,label=L"\textbf{M}\;",legendfontsize=15,legend=:outertop,legend_column=-1,
tickfontsize=12,xlabel=L"Species\; response \;correlation \;(q)",foreground_color_legend=nothing,
ylab=L"Fluctuation\; size",guidefontsize=18,yaxis=(formatter=y->string(round( y / 10^-4;digits=3))),c=:blue,
y_foreground_color_axis=:blue,y_foreground_color_border=:blue,xlim=(0,1.0))
annotate!([(0.0,maximum(pp[:,1])*1.008,text(L"\times10^{-4}", 14, :black))])
annotate!([(1.04,maximum(pp[:,1])*1.008,text(L"\times10^{-2}", 14, :black))])

plot!(twinx(),q,pp[:,2],linewidth=5,color=:red,widen=false,grid=false,xlim=(0,1.0),
yaxis=(formatter=y->string(round( y / 10^-2;digits=3))),legendfontsize=11,legend=:none,
tickfontsize=12,y_foreground_color_axis=:red,y_foreground_color_border=:red,label="",dpi=750)

plot!(twiny(),xticks=:false,yticks=:false,y_foreground_color_border=:blue)
plot!(legend=(0.45,1.05),top_margin=10mm)

savefig("coral_std_m3_q.png")


####################################------Plot flactuation size w.r.t tempral environmental correlation------###################


pp=readdlm("coral_sd_trend_k.txt")
q=0
k=[0:0.1:0.9;]

plot([pp[1,1]],label=L"\textbf{C}\;\;\;",c=:red,lw=5,foreground_color_legend=nothing,legendfontsize=15,legend=:outertop)
plot!(k,pp[:,1],linewidth=5,box=:on,widen=false,grid=false,label=L"\textbf{M}\;",legendfontsize=15,legend=:outertop,legend_column=-1,
tickfontsize=12,xlabel=L"Temporal \;correlation \;(k)",foreground_color_legend=nothing,
ylab=L"Fluctuation\; size",guidefontsize=18,yaxis=(formatter=y->string(round( y / 10^-3;digits=3))),c=:blue,
y_foreground_color_axis=:blue,y_foreground_color_border=:blue,xlim=(0,0.9))
annotate!([(0.0,maximum(pp[:,1])*1.05,text(L"\times10^{-3}", 14, :black))])
annotate!([(0.93,maximum(pp[:,1])*1.05,text(L"\times10^{-2}", 14, :black))])

plot!(twinx(),k,pp[:,2],linewidth=5,color=:red,widen=false,grid=false,xlim=(0,0.9),
yaxis=(formatter=y->string(round( y / 10^-2;digits=3))),legendfontsize=11,legend=:none,
tickfontsize=12,y_foreground_color_axis=:red,y_foreground_color_border=:red,label="",dpi=750)

plot!(twiny(),xticks=:false,yticks=:false,y_foreground_color_border=:blue)
plot!(legend=(0.45,1.05),top_margin=10mm)

savefig("coral_std_m3_k.png")







###############################----------Plot flactuation size for general model of two species mutualism----------------######################### 


####################################------Plot flactuation size w.r.t species response correlation------###################

q=[0:0.1:1;]
k=0
pp=readdlm("mod4_ob_sd_trend_q.txt")

plot([pp[1,1]],label=L"\textbf{M_2}\;\;\;",c=:red,lw=5,foreground_color_legend=nothing,legendfontsize=15,legend=:outertop)
plot!(q,pp[:,1],linewidth=8,box=:on,widen=false,grid=false,label=L"\textbf{M_1}\;",legendfontsize=15,legend=:outertop,legend_column=-1,
tickfontsize=12,xlabel=L"Species\; response \;correlation \;(q)",foreground_color_legend=nothing,
ylab=L"Fluctuation\; size",guidefontsize=18,yaxis=(formatter=y->string(round( y / 10^-2;digits=3))),c=:blue,
y_foreground_color_axis=:blue,y_foreground_color_border=:blue,xlim=(0,1.0))
annotate!([(0.0,maximum(pp[:,1])*1.01,text(L"\times10^{-2}", 14, :black))])
annotate!([(1.04,maximum(pp[:,1])*1.01,text(L"\times10^{-2}", 14, :black))])

plot!(twinx(),q,pp[:,2],linewidth=5,color=:red,widen=false,grid=false,xlim=(0,1.0),
yaxis=(formatter=y->string(round( y / 10^-2;digits=3))),legendfontsize=11,legend=:none,
tickfontsize=12,y_foreground_color_axis=:red,y_foreground_color_border=:red,label="",dpi=750)

plot!(twiny(),xticks=:false,yticks=:false,y_foreground_color_border=:blue)
plot!(legend=(0.45,1.05),top_margin=10mm)

savefig("trend_std_m4_ob_q.png")


####################################------Plot flactuation size w.r.t tempral environmental correlation------###################


pp=readdlm("mod4_ob_sd_trend_k.txt")
q=0
k=[0:0.1:0.9;]

plot([pp[1,1]],label=L"\textbf{M_2}\;\;\;",c=:red,lw=5,foreground_color_legend=nothing,legendfontsize=15,legend=:outertop)
plot!(k,pp[:,1],linewidth=8,box=:on,widen=false,grid=false,label=L"\textbf{M_1}\;",legendfontsize=15,legend=:outertop,legend_column=-1,
tickfontsize=12,xlabel=L"Temporal \;correlation \;(k)",foreground_color_legend=nothing,
ylab=L"Fluctuation\; size",guidefontsize=18,yaxis=(formatter=y->string(round( y / 10^-2;digits=3))),c=:blue,
y_foreground_color_axis=:blue,y_foreground_color_border=:blue,xlim=(0,0.9))
annotate!([(0.0,maximum(pp[:,1])*1.05,text(L"\times10^{-2}", 14, :black))])
annotate!([(0.93,maximum(pp[:,1])*1.05,text(L"\times10^{-2}", 14, :black))])

plot!(twinx(),k,pp[:,2],linewidth=5,color=:red,widen=false,grid=false,xlim=(0,0.9),
yaxis=(formatter=y->string(round( y / 10^-2;digits=3))),legendfontsize=11,legend=:none,
tickfontsize=12,y_foreground_color_axis=:red,y_foreground_color_border=:red,label="",dpi=750)

plot!(twiny(),xticks=:false,yticks=:false,y_foreground_color_border=:blue)
plot!(legend=(0.45,1.05),top_margin=10mm)

savefig("trend_std_m4_ob_k.png")





###############################----------Plot flactuation size for Plant-Pollinator model----------------######################### 


####################################------Plot flactuation size w.r.t species response correlation------###################


q=[0:0.1:1;]
k=0
pp=readdlm("mod5_bf_sd_trend_q.txt")

plot([pp[1,1]],label=L"\textbf{A}\;\;\;",c=:red,lw=5,foreground_color_legend=nothing,legendfontsize=15,legend=:outertop)
plot!(q,pp[:,1],linewidth=8,box=:on,widen=false,grid=false,label=L"\textbf{P}\;",legendfontsize=15,legend=:outertop,legend_column=-1,
tickfontsize=12,xlabel=L"Species\; response \;correlation \;(q)",foreground_color_legend=nothing,
ylab=L"Fluctuation\; size",guidefontsize=18,yaxis=(formatter=y->string(round( y / 10^-3;digits=3))),c=:blue,
y_foreground_color_axis=:blue,y_foreground_color_border=:blue,xlim=(0,1.0))
annotate!([(0.0,maximum(pp[:,1])*1.005,text(L"\times10^{-3}", 14, :black))])
annotate!([(1.04,maximum(pp[:,1])*1.005,text(L"\times10^{-3}", 14, :black))])

plot!(twinx(),q,pp[:,2],linewidth=5,color=:red,widen=false,grid=false,xlim=(0,1.0),
yaxis=(formatter=y->string(round( y / 10^-3;digits=3))),legendfontsize=11,legend=:none,
tickfontsize=12,y_foreground_color_axis=:red,y_foreground_color_border=:red,label="",dpi=750)

plot!(twiny(),xticks=:false,yticks=:false,y_foreground_color_border=:blue)
plot!(legend=(0.45,1.05),top_margin=10mm)

savefig("trend_std_m5_fac_q.png")


####################################------Plot flactuation size w.r.t tempral environmental correlation------###################


pp=readdlm("mod5_bf_sd_trend_k.txt")
q=0
k=[0:0.1:0.9;]

plot([pp[1,1]],label=L"\textbf{A}\;\;\;",c=:red,lw=5,foreground_color_legend=nothing,legendfontsize=15,legend=:outertop)
plot!(k,pp[:,1],linewidth=8,box=:on,widen=false,grid=false,label=L"\textbf{P}\;",legendfontsize=15,legend=:outertop,legend_column=-1,
tickfontsize=12,xlabel=L"Temporal \;correlation \;(k)",foreground_color_legend=nothing,
ylab=L"Fluctuation\; size",guidefontsize=18,yaxis=(formatter=y->string(round( y / 10^-1;digits=3))),c=:blue,
y_foreground_color_axis=:blue,y_foreground_color_border=:blue,xlim=(0,0.9))
annotate!([(0.0,maximum(pp[:,1])*1.07,text(L"\times10^{-1}", 14, :black))])
annotate!([(0.93,maximum(pp[:,1])*1.07,text(L"\times10^{-1}", 14, :black))])

plot!(twinx(),k,pp[:,2],linewidth=5,color=:red,widen=false,grid=false,xlim=(0,0.9),
yaxis=(formatter=y->string(round( y / 10^-2;digits=3))),legendfontsize=11,legend=:none,
tickfontsize=12,y_foreground_color_axis=:red,y_foreground_color_border=:red,label="",dpi=750)

plot!(twiny(),xticks=:false,yticks=:false,y_foreground_color_border=:blue)
plot!(legend=(0.45,1.05),top_margin=10mm)

savefig("trend_std_m5_fac_k.png")



###############################----------Plot flactuation size for patch occupancy model for mutualism----------------######################### 


####################################------Plot flactuation size w.r.t species response correlation------###################

q=[0:0.1:1;]
k=0
pp=readdlm("mod6_spmut_sd_trend_q.txt")

plot([pp[1,1]],label=L"\textbf{a}\;\;\;",c=:red,lw=5,foreground_color_legend=nothing,legendfontsize=15,legend=:outertop)
plot!(q,pp[:,1],linewidth=5,box=:on,widen=false,grid=false,label=L"\textbf{p}\;",legendfontsize=15,legend=:outertop,legend_column=-1,
tickfontsize=12,xlabel=L"Species\; response \;correlation \;(q)",foreground_color_legend=nothing,
ylab=L"Fluctuation\; size",guidefontsize=18,yaxis=(formatter=y->string(round( y / 10^-2;digits=3))),c=:blue,
y_foreground_color_axis=:blue,y_foreground_color_border=:blue,xlim=(0,1.0))
annotate!([(0.0,maximum(pp[:,1])*1.003,text(L"\times10^{-2}", 14, :black))])
annotate!([(1.04,maximum(pp[:,1])*1.003,text(L"\times10^{-2}", 14, :black))])

plot!(twinx(),q,pp[:,2],linewidth=5,color=:red,widen=false,grid=false,xlim=(0,1.0),
yaxis=(formatter=y->string(round( y / 10^-2;digits=3))),legendfontsize=11,legend=:none,
tickfontsize=12,y_foreground_color_axis=:red,y_foreground_color_border=:red,label="",dpi=750)

plot!(twiny(),xticks=:false,yticks=:false,y_foreground_color_border=:blue)
plot!(legend=(0.45,1.05),top_margin=10mm)

savefig("trend_std_m6_spmut_q.png")


####################################------Plot flactuation size w.r.t tempral environmental correlation------###################


pp=readdlm("mod6_spmut_sd_trend_k.txt")
q=0
k=[0:0.1:0.9;]

plot([pp[1,1]],label=L"\textbf{a}\;\;\;",c=:red,lw=5,foreground_color_legend=nothing,legendfontsize=15,legend=:outertop)
plot!(k,pp[:,1],linewidth=5,box=:on,widen=false,grid=false,label=L"\textbf{p}\;",legendfontsize=15,legend=:outertop,legend_column=-1,
tickfontsize=12,xlabel=L"Temporal \;correlation \;(k)",foreground_color_legend=nothing,
ylab=L"Fluctuation\; size",guidefontsize=18,yaxis=(formatter=y->string(round( y / 10^-2;digits=3))),c=:blue,
y_foreground_color_axis=:blue,y_foreground_color_border=:blue,xlim=(0,0.9))
annotate!([(0.0,maximum(pp[:,1])*1.05,text(L"\times10^{-2}", 14, :black))])
annotate!([(0.93,maximum(pp[:,1])*1.05,text(L"\times10^{-2}", 14, :black))])

plot!(twinx(),k,pp[:,2],linewidth=5,color=:red,widen=false,grid=false,xlim=(0,0.9),
yaxis=(formatter=y->string(round( y / 10^-2;digits=3))),legendfontsize=11,legend=:none,
tickfontsize=12,y_foreground_color_axis=:red,y_foreground_color_border=:red,label="",dpi=750)

plot!(twiny(),xticks=:false,yticks=:false,y_foreground_color_border=:blue)
plot!(legend=(0.45,1.05),top_margin=10mm)

savefig("trend_std_m6_spmut_k.png")


