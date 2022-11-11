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
set_default(fmt = "%.5f",           # latexify
    convert_unicode = false,
    latex = false)


#======================================================#
#       TABLE + FIGURE
#------------------------------------------------------#
# Produce Table 1
function estimates_table(; output1::Vector{Float64}, output2::Vector{Float64}, cf::String)
     # define the latex table head and foot
    tabhead = "
    \\begin{table}\\caption{Estimates: efficient and inefficient }\n\\centering
    \\begin{tabular}{|l|c|c|}
    \\toprule
    \tParameter & b1 & b2 \\\\
    \\miderule"

      tabfoot = "
      \\bottomrule
      \\end{tabular}
      \\end{table}"

    # use latexify to create the table part (minus the column names, those are too complicated for latexify)
    tabside   = [ "\trho", "\tsigma"]
    tabvals   = hcat(output1, output2)
    temptable = latexify(tabvals, env = :tabular, side = tabside)

    # get rid of latexify's table head and foot, so that we can replace it with our own.
    temptable = temptable[25:end-14]
    # put everything together and write the table
    write(joinpath(figpath, "est_"*cf*"_.tex"), tabhead*temptable*tabfoot)
end

function errors_table(; output1::Vector{Float64}, output2::Vector{Float64}, cf::String)
     # define the latex table head and foot
    tabhead = "
    \\begin{table}\\caption{Standard errors}\n\\centering
    \\begin{tabular}{|l|c|c|}
    \\toprule
    \tParameter & Estimate & St. error \\\\
    \\miderule"

      tabfoot = "
      \\bottomrule
      \\end{tabular}
      \\end{table}"

    # use latexify to create the table part (minus the column names, those are too complicated for latexify)
    tabside   = [ "\trho", "\tsigma"]
    tabvals   = hcat(output1, output2)
    temptable = latexify(tabvals, env = :tabular, side = tabside)

    # get rid of latexify's table head and foot, so that we can replace it with our own.
    temptable = temptable[25:end-14]
    # put everything together and write the table
    write(joinpath(figpath, "err_"*cf*"_.tex"), tabhead*temptable*tabfoot)
end

function jacobian_table(; output1::Matrix{Float64}, cf::String)
     # define the latex table head and foot
    tabhead = "
    \\begin{table}\\caption{Jacobian approximation}\n\\centering
    \\begin{tabular}{|l|c|c|}
    \\toprule
    \\miderule"

      tabfoot = "
      \\bottomrule
      \\end{tabular}
      \\end{table}"

    # use latexify to create the table part (minus the column names, those are too complicated for latexify)
    temptable = latexify(output1, env = :tabular)

    # get rid of latexify's table head and foot, so that we can replace it with our own.
    temptable = temptable[25:end-14]
    # put everything together and write the table
    write(joinpath(figpath, "jacobian_"*cf*"_.tex"), tabhead*temptable*tabfoot)
end

function j_test_table(; output1::Float64, output2::Float64, output3::Float64)
     # define the latex table head and foot
    tabhead = "
    \\begin{table}\\caption{J test}\n\\centering
    \\begin{tabular}{|l|c|c|c|}
    \\toprule
    \t J test statistic & mean + var & var + autocor & mean + var + autocor  \\\\
    \\miderule"

      tabfoot = "
      \\bottomrule
      \\end{tabular}
      \\end{table}"

    # use latexify to create the table part (minus the column names, those are too complicated for latexify)
    temptable = latexify(hcat(output1, output2, output3), env = :tabular)

    # get rid of latexify's table head and foot, so that we can replace it with our own.
    temptable = temptable[25:end-14]
    # put everything together and write the table
    write(joinpath(figpath, "j_test_res.tex"), tabhead*temptable*tabfoot)
end
