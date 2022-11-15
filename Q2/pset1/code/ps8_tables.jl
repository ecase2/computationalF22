#=
    * PROJECT:     COMPUTATIONAL FALL 2022 PSET 6
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
function table_score(; output1::Array{Float64}, output2::Array{Float64})
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

function table_hessian_diff1(; difference::Matrix{Float64}, name::String)
     # define the latex table head and foot
    tabhead = "
    \\begin{table}\\caption{Absolute difference between the score of the log-likelihood function and numerical derivative}\n\\centering
    \\begin{sideways}
    \\begin{tabular}{|c|c|c|c|c|c|c|c|}
    \\toprule
    \t 1 & 2& 3 & 4 & 5 & 6 & 7 & 8 \\\\
    \\miderule"

      tabfoot = "
      \\bottomrule
      \\end{tabular}
      \\end{sideways}
      \\end{table}"

    # use latexify to create the table part (minus the column names, those are too complicated for latexify)
    tabvals   = difference
    temptable = latexify(tabvals, env = :tabular; fmt="%.1e")

    # get rid of latexify's table head and foot, so that we can replace it with our own.
  #  temptable = temptable[25:end-14]

    # put everything together and write the table
    write(joinpath(tabpath, name), tabhead*temptable*tabfoot)
end

function table_hessian_diff2(; difference::Matrix{Float64}, name::String)
     # define the latex table head and foot
    tabhead = "
    \\begin{table}\\caption{Absolute difference between the score of the log-likelihood function and numerical derivative}\n\\centering
    \\begin{sideways}
    \\begin{tabular}{|c|c|c|c|c|c|c|c|c|}
    \\toprule
    \t 9 & 10 & 11 & 12 & 13 & 14 & 15 & 16 & 17 \\\\
    \\miderule"

      tabfoot = "
      \\bottomrule
      \\end{tabular}
      \\end{sideways}
      \\end{table}"

    # use latexify to create the table part (minus the column names, those are too complicated for latexify)
    tabvals   = difference
    temptable = latexify(tabvals, env = :tabular; fmt="%.1e")

    # get rid of latexify's table head and foot, so that we can replace it with our own.
  #  temptable = temptable[25:end-14]

    # put everything together and write the table
    write(joinpath(tabpath, name), tabhead*temptable*tabfoot)
end

function estimates_Newton(; output1::Vector{Float64})
    tabhead = "
    \\begin{table}\\caption{Estimates obtained by the Newton algorithm}\n\\centering
    \\begin{tabular}{|l|c|}
    \\toprule
    \t Regressor & Newton \\\\
    \\miderule"

      tabfoot = "
      \\bottomrule
      \\end{tabular}
      \\end{table}"
    # use latexify to create the table part (minus the column names, those are too complicated for latexify)
    tabside   = [ "\ti_large_loan", "\ti_medium_loan", "\trate_spread", "\ti_refinance", "\tage_r", "\tcltv", "\tdti", "\tcu", "\tfirst_mort_r",
    "\tscore_0", "\tscore_1", "\ti_FHA",  "\ti_open_year2", "\ti_open_year3", "\ti_open_year4", "\ti_open_year5", "\tconstant"]

    tabvals   = output1
    temptable = latexify(tabvals, env = :tabular)

    # get rid of latexify's table head and foot, so that we can replace it with our own.
  #  temptable = temptable[25:end-14]

    # put everything together and write the table
    write(joinpath(tabpath, "ps8_results_Newton.tex"), tabhead*temptable*tabfoot)
end

function estimates(; output1::Array{Float64}, output2::Array{Float64}, output3::Array{Float64})
    tabhead = "
    \\begin{table}\\caption{Comparison of the estimates}\n\\centering
    \\begin{tabular}{|l|c|c|c|}
    \\toprule
    \t Regressor & Newton & BFGS & Simplex \\\\
    \\miderule"

      tabfoot = "
      \\bottomrule
      \\end{tabular}
      \\end{table}"
    # use latexify to create the table part (minus the column names, those are too complicated for latexify)
    tabside   = [ "\ti_large_loan", "\ti_medium_loan", "\trate_spread", "\ti_refinance", "\tage_r", "\tcltv", "\tdti", "\tcu", "\tfirst_mort_r",
    "\tscore_0", "\tscore_1", "\ti_FHA", "\ti_open_year1", "\ti_open_year2", "\ti_open_year3", "\ti_open_year4", "\ti_open_year5", "\tconstant"]

    tabvals   = hcat(output1, output2, output3)
    temptable = latexify(tabvals, env = :tabular; side = tabside)

    # get rid of latexify's table head and foot, so that we can replace it with our own.
  #  temptable = temptable[25:end-14]

    # put everything together and write the table
    write(joinpath(tabpath, "ps8_results.tex"), tabhead*temptable*tabfoot)
end
