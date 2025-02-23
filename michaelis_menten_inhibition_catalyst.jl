using Catalyst
using DifferentialEquations
using Plots

# Define the reaction network
rn = @reaction_network begin
    k1, S + E ↔ ES, k2
    k3, ES → E + P
    k4, I + E ↔ EI, k5
end k1 k2 k3 k4 k5

# Initial concentrations
u0 = [S => 10.0, E => 1.0, ES => 0.0, P => 0.0, I => 1.0, EI => 0.0]

# Parameters: k1, k2, k3, k4, k5
p = [0.1, 0.1, 0.1, 0.1, 0.1]

# Time span
tspan = (0.0, 50.0)

# Solve the differential equations
prob = ODEProblem(rn, u0, tspan, p)
sol = solve(prob)

# Plot the results
plot(sol, xlabel="Time", ylabel="Concentration", label=["S" "E" "ES" "P" "I" "EI"])