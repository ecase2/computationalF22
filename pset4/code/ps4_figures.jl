#=
    * PROJECT:     COMPUTATIONAL FALL 2022 PSET 4
    * AUTHORS:     Hanna Han, Emily Case, Anna Lukianova
    * CONTENTS:    Creates functions for write up graphs/tables 

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
using Plots

# set defaults:
default(titlefont = (20, "times"),  # plots
    xlabel    = "time", 
    linewidth = 2)
#======================================================#



#======================================================#
#       TRANSITION PATH GRAPHS
#------------------------------------------------------#

function graphPath(par::parameters, res::results)
    @unpack T = par 
    x = 1:T
    #! need to save exercise 2 graphs differently. 
    # interest rate 
    rPath = plot(x, res.r, title = "Interest rate transition path", ylabel = "r")
    savefig(joinpath(figpath, "rPath.png"))

    # wage rate 
    wPath = plot(x, res.w, title = "Wages transition path", ylabel = "w")
    savefig(joinpath(figpath, "wPath.png"))

    # aggregate labor 
    LPath = plot(x, res.L, title = "Labor transition path", ylabel = "L")
    savefig(joinpath(figpath, "LPath.png"))

    # aggregate capital 
    KPath = plot(x, res.K, title = "Capital transition path", ylabel = "K")
    savefig(joinpath(figpath, "KPath.png"))

    return rPath, wPath, LPath, KPath
end


#======================================================#
