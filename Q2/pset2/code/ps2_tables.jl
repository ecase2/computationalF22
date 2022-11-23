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
function table_comparison(; output1::Float64, output2::Float64, output3::Float64, name::String)
     # define the latex table head and foot
    tabhead = "
    \\begin{table}\\caption{Predcited values of the log-likelihood computed using 3 approaches}\n\\centering
    \\begin{tabular}{|c|c|c|}
    \\toprule
    \t Quadrature & GHK & Accept/Reject\\\\
    \\miderule"

      tabfoot = "
      \\bottomrule
      \\end{tabular}
      \\end{table}"
    # use latexify to create the table part (minus the column names, those are too complicated for latexify)
    tabvals   = hcat(output1, output2, output3)
    temptable = latexify(tabvals, env = :tabular)

    # put everything together and write the table
    write(joinpath(tabpath, name), tabhead*temptable*tabfoot)
end
