#=
    * THIS FILE WRITES FUNCTIONS FOR ALL OF OUR GRAPHS & TABLES 

    * NOTE: this file only writes *functions* which, when called, will create and save figures.
    *       the functions should be called in other files, in order to produce figures. the reasoning
    *       for this method is (1) that julia is faster when things are in functions, and (2) gives 
    *       us more flexibility about when in our program we need to create graphs. 

    * NOTE: when saving graphs, use savefig(joinpath(figpath, "nameofgraph.png"))
    *       when saving tex files, use write(joinpath(figpath, "nameoffile.tex"), stringobject)
    *       
    *       using joinpath() and figpath allows us to always save figures to the same directory,
    *       regradless of who's computer it is on (mac or windows). you must have figpath defined
    *       in order for it to work... check that figpath routes to pset3/code/figs! if it doesn't 
    *       look right, your working directory might be off. check that pwd() returns the file path
    *       to computationalF22 (the whole repository)
=#

#======================================================#
#       SET UP: IMPORT PACKAGES, SET DEFAULTS 
#------------------------------------------------------#
using Plots, Latexify

# set defaults:
default(titlefont = (20, "times"),  # plots
    linewidth = 2)
set_default(fmt = "%.3f",           # latexify 
    convert_unicode = false, 
    latex = false)  
#======================================================#



#======================================================#
#       EXERCISE 1 GRAPHS
#------------------------------------------------------#
function createAllGraphs(par::parameters, res::results)
    @unpack pol_func, val_func = res
    @unpack a_grid = par

    # savings of the worker of 20 years old with high ans low productivity shocks
    # Anna: I don't like the policy function
    plot(a_grid, pol_func[:, 1, 20],         
        label = "High productivity")
    plot!(a_grid, pol_func[:, 2, 20], 
        label = "Low productivity")

    # labels, title, legend, etc.: 
    plot!(xlabel = "Current assets, a", 
        ylabel   = "Future assets, a", 
        title    = "Savings of the worker at 20 years old", 
        legend   = :bottomright)

    savefig(joinpath(figpath, "pol_func.png"))

    # # value function of a retired of 50 years old
    # plot(a_grid, val_func[:, 1, 50], xlabel = "Current assets, a", ylabel = "Value", legend = false, title = "Value function of the retired at the age 50")
    # savefig(joinpath(figpath, "savings_at_20.png"))

    # value function of a retired of 50 years old
    plot(a_grid, val_func[:, 1, 50], xlabel = "Current assets, a", ylabel = "Value", title = "Value function of the retired at the age 50")
    savefig(joinpath(figpath, "value_at_50.png"))

    # plot(a_grid, pol_func[:, 1, 20], label = "High productivity", xlabel = "Current assets, a", ylabel = "Future assets, a", title = "Savings of the worker at 20 years old")
    # plot!(a_grid, a_grid, label = "45 degree line", linestyle = :dash)
    # plot!(a_grid, pol_func[:, 2, 20], label = "Low productivity", legend =:bottomright)

    # savefig(joinpath(figpath, "savings_at_20.png"))

    # value function of a retired of 50 years old
    plot(a_grid, val_func[:, 1, 50], xlabel = "Current assets, a", ylabel = "Value", title = "Value function of the retired at the age 50")
    savefig(joinpath(figpath,"value_at_50.png"))
end
#======================================================#



#======================================================#
#       EXERCISE 3 TABLE  
#------------------------------------------------------#
function resultsTable(bm_ss, bm_noss, noshock_ss, noshock_noss, inelasticl_ss, inelasticl_noss)
    # define the latex table head and foot 
    tabhead = "
    \\begin{table}\\caption{Exercise 3 Results}\n\\centering
    \\begin{tabular}{lcccccc}
    \\toprule
    \t& \\multicolumn{2}{c}{Benchmark model} & \\multicolumn{2}{c}{No risk, \$z^L = z^H = 0.5\$} & \\multicolumn{2}{c}{Exogenous labor, \$\\gamma = 1\$} \\\\ 
    \t\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}
    \t& with SS & wo SS & with SS & wo SS & with SS & wo SS\\\\
    \\midrule"

    tabfoot = "
    \\bottomrule
    \\end{tabular}
    \\end{table}"

    # use latexify to create the table part (minus the column names, those are too complicated for latexify) 
    tabside   = [ "\tcapital \$K\$", "\tlabor \$L\$", "\twage \$w\$", "\tinterest \$r\$", "\tpension benefit \$b\$", "\ttotal welfare \$W\$", "\tcv(wealth)"]
    tabvals   = hcat(bm_ss, bm_noss, noshock_ss, noshock_noss, inelasticl_ss, inelasticl_noss)
    temptable = latexify(tabvals, env = :tabular, side = tabside)

    # get rid of latexify's table head and foot, so that we can replace it with our own. 
    temptable = temptable[25:end-14]
    # put everything together and write the table 
    write(joinpath(figpath, "resultstable.tex"), tabhead*temptable*tabfoot)
end
