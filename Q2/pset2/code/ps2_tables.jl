#=
    * PROJECT:     COMPUTATIONAL FALL 2022 PSET 2 (Q2)
    * AUTHORS:     Hanna Han and Anna Lukianova
    * CONTENTS:    Creates functions for write up graphs/tables
=#

#======================================================#
#       SET UP: IMPORT PACKAGES, SET DEFAULTS
#------------------------------------------------------#
# set defaults:
#default(titlefont = (20, "times"),  # plots
#    linewidth = 2)
set_default(fmt = "%.3f",           # latexify
    convert_unicode = false,
    latex = false)


#======================================================#
#       TABLE + FIGURE
#------------------------------------------------------#
# Produce Table 1
# NEED TO WRITE THE FUNCTION TO CREATE A TABLE FOR TASK 4 (comparison table)
function table_comparison(; output1::Array{Float64}, output2::Array{Float64})
     # define the latex table head and foot
    tabhead = "
    \\begin{table}\\caption{Comparison of the score of the log-likelihood function and numerical derivative}\n\\centering
    \\begin{tabular}{|c|c|}
    \\toprule
    \t Score & Derivative\\\\
    \\miderule"

      tabfoot = "
      \\bottomrule
      \\end{tabular}
      \\end{table}"
    # use latexify to create the table part (minus the column names, those are too complicated for latexify)
    tabvals   = hcat(output1, output2)
    temptable = latexify(tabvals, env = :tabular)

    # get rid of latexify's table head and foot, so that we can replace it with our own.
  #  temptable = temptable[25:end-14]

    # put everything together and write the table
    write(joinpath(tabpath, "ps8_results_score.tex"), tabhead*temptable*tabfoot)
end
