#=
    * PROJECT:     COMPUTATIONAL FALL 2022 PSET 6
    * AUTHORS:     Hanna Han and Anna Lukianova
    * CONTENTS:    Creates functions for write up graphs/tables
=#

#======================================================#
#       SET UP: IMPORT PACKAGES, SET DEFAULTS
#------------------------------------------------------#
# set defaults:
default(titlefont = (20, "times"),  # plots
    linewidth = 2)
set_default(fmt = "%.3f",           # latexify 
    convert_unicode = false, 
    latex = false)  


#======================================================#
#       TABLE + FIGURE
#------------------------------------------------------#
# Produce Table 1
function write_table(output1, output2, output3; cf = "")
     # define the latex table head and foot 
    tabhead = "
    \\begin{table}\\caption{Final model moments}\n\\centering
    \\begin{tabular}{lccc}
    \\toprule
    \t Moments & Standard & TV1 shock (\$\alpha=1\$) & TV1 shock (\$\alpha=2\$) \\\\
    \\miderule"

    tabfoot = "
    \\bottomrule
    \\end{tabular}
    \\end{table}"

    # use latexify to create the table part (minus the column names, those are too complicated for latexify) 
    tabside   = [ "\tPrice level", "\tMass of incumbents", "\tMass of entrants", "\tMass of exits", "\tAggregate labor", "\tLabor of incumbents", "\tLabor of entrants", "\tFraction of labor in entrants"]
    tabvals   = hcat(output1, output2, output3)
    temptable = latexify(tabvals, env = :tabular, side = tabside)

    # get rid of latexify's table head and foot, so that we can replace it with our own. 
    temptable = temptable[25:end-14]
    # put everything together and write the table 
    write(joinpath(figpath, "ps6_results"*cf*".tex"), tabhead*temptable*tabfoot)
end

# Produce Figure 1
function graphExit(par::Params, pol_func1, pol_func2, pol_func3)
    @unpack Ns = par

    x = 1:Ns

    plot(x, pol_func1, title = "Probability of market exit by firm type", xlabel = "Firm type", ylabel = "Probability of exit", label = "Standard")
    plot!(x, pol_func2, label = "TV1 shock \$\\alpha = 1\$")
    plot!(x, pol_func3, label = "TV1 shock \$\\alpha = 2\$")

    savefig(joinpath(figpath, "fig1.png"))
end