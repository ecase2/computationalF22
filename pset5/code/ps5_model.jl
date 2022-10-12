#=
    PROJECT:     COMPUTATIONAL FALL 2022 PSET 5
    AUTHORS:     Hanna Han, Emily Case, Anna Lukianova
    CONTENTS:    This file contains main functions.
=#

# Initialize
function Initialize()
    par = Params()
end

# Simulate aggregate capital path
function simulate_KPath()

end

# Estimate regression
function estimate_reg()
    # Sort based on z

    # Estimate regression

end

# Solve model
function solve_model(; λ = 0.5, iter = 1000, tol = 0.005)
    # Draw shocks

    # Set initial guess
    a0_old = 0.095
    a1_old = 0.999
    b0_old = 0.085
    b1_old = 0.999

    diff = 10.0
    n = 1

    # 
    while (diff > tol && n < iter)
        println("BEGINNING ITERATION $n")
        n = n+1

        # Solve

        # Simulate capital path

        # Estimate regression
        a0_new, a1_new, b0_new, b1_new = estimate_reg()

        # Calculate difference
        diff = abs(a0_new - a0_old) + abs(a1_new - a1_old) + abs(b0_new - b0_old) + abs(b1_new - b1_old)

        # Adjust guess
        if diff > tol
            a0_old = λ*a0_new + (1-λ)*a0_old
            a1_old = λ*a1_new + (1-λ)*a1_old
            b0_old = λ*b0_new + (1-λ)*b0_old
            b1_old = λ*b1_new + (1-λ)*b1_old
        else
            #[Insert something here]

            println("DONE!\n")
            break
        end
    end

    return 
end