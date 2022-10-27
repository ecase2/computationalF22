#=
    * PROJECT:     COMPUTATIONAL FALL 2022 PSET 6
    * AUTHORS:     Hanna Han and Anna Lukianova
    * CONTENTS:    Creates functions for write up graphs/tables
=#

#======================================================#
#       SET UP: IMPORT PACKAGES, SET DEFAULTS
#------------------------------------------------------#
using Plots

# set defaults:
default(titlefont = (20, "times"),  # plots
    xlabel    = "time",
    linewidth = 2)
set_default(fmt = "%.3f",           # latexify 
    convert_unicode = false, 
    latex = false) 


#======================================================#
#       TABLE + FIGURE
#------------------------------------------------------#
# Produce Table 1
function write_table(output)
     # define the latex table head and foot 
    tabhead = "
    \\begin{table}\\caption{Final model moments}\n\\centering
    \\begin{tabular}{lccc}
    \\toprule
    \t Variable & Standard & TV1 shock (\$\alpha=1\$) & TV1 shock (\$\alpha=2\$) \\\\
    \\miderule"

    tabfoot = "
    \\bottomrule
    \\end{tabular}
    \\end{table}"

    # use latexify to create the table part (minus the column names, those are too complicated for latexify) 
    tabside   = [ "\tPrice level", "\tMass of incumbents", "\tMass of entrants", "\tMass of exits", "\tAggregate labor", "\tLabor of incumbents", "\tLabor of entrants", "\tFraction of labor in entrants"]
    tabvals   = output
    temptable = latexify(tabvals, env = :tabular, side = tabside)

    # get rid of latexify's table head and foot, so that we can replace it with our own. 
    temptable = temptable[25:end-14]
    # put everything together and write the table 
    write(joinpath(figpath, "ps6_results.tex"), tabhead*temptable*tabfoot)
end

# Produce Figure 1
function graphPath(prim::primitives, res::TransitionResults; cf="")
    @unpack T = prim
    x = 1:T

    # interest rate
    rPath = plot(x, res.r, title = "Interest rate transition path", ylabel = "r")
    savefig(joinpath(figpath, "rPath"*cf*".png"))

    # wage rate
    wPath = plot(x, res.w, title = "Wages transition path", ylabel = "w")
    savefig(joinpath(figpath, "wPath"*cf*".png"))

    # aggregate labor
    LPath = plot(x, res.L, title = "Labor transition path", ylabel = "L")
    savefig(joinpath(figpath, "LPath"*cf*".png"))

    # aggregate capital
    KPath = plot(x, res.K, title = "Capital transition path", ylabel = "K")
    savefig(joinpath(figpath, "KPath"*cf*".png"))

    return rPath, wPath, LPath, KPath
end