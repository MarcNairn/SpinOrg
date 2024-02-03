using LsqFit

# a two-parameter exponential model
# x: array of independent variables
# p: array of model parameters
# model(x, p) will accept the full data set as the first argument `x`.
# This means that we need to write our model function so it applies
# the model to the full dataset. We use `@.` to apply the calculations
# across all rows.
# @. model(x, p) = p[1]*exp(-x*p[2])


# # # Define the q-Gaussian model
# # @. q_exponential(x, q) = (1 + (1 - q) * x)^(1 / (1 - q))
# # @. q_gaussian(x, p) = p[1] * q_exponential(-p[2] * x^2, p[3])

# # some example data
# # xdata: independent variables
# # ydata: dependent variable
# xdata = range(0, stop=10, length=20)
# ydata = model(xdata, [1.0 2.0]) + 0.01*randn(length(xdata))
# p0 = [0.5, 0.5]

# fit = curve_fit(model, xdata, ydata, p0)
# # fit is a composite type (LsqFitResult), with some interesting values:
# #	dof(fit): degrees of freedom
# #	coef(fit): best fit parameters
# #	fit.resid: residuals = vector of residuals
# #	fit.jacobian: estimated Jacobian at solution
# lb = [1.1, -0.5]
# ub = [1.9, Inf]
# p0_bounds = [1.2, 1.2] # we have to start inside the bounds
# # Optional upper and/or lower bounds on the free parameters can be passed as an argument.
# # Bounded and unbouded variables can be mixed by setting `-Inf` if no lower bounds
# # is to be enforced for that variable and similarly for `+Inf`
# fit_bounds = curve_fit(model, xdata, ydata, p0_bounds, lower=lb, upper=ub)

# # We can estimate errors on the fit parameters,
# # to get standard error of each parameter:
# sigma = stderror(fit)
# # to get margin of error and confidence interval of each parameter at 5% significance level:
# margin_of_error = margin_error(fit, 0.05)
# confidence_inter = confint(fit; level=0.95)

# # The finite difference method is used above to approximate the Jacobian.
# # Alternatively, a function which calculates it exactly can be supplied instead.
# function jacobian_model(x,p)
#     J = Array{Float64}(undef, length(x), length(p))
#     @. J[:,1] = exp(-x*p[2])     #dmodel/dp[1]
#     @. @views J[:,2] = -x*p[1]*J[:,1] #dmodel/dp[2], thanks to @views we don't allocate memory for the J[:,1] slice
#     J
# end
# fit = curve_fit(model, jacobian_model, xdata, ydata, p0)


# Define the q-Gaussian model
@. q_exponential(x, q) = (1 + (1 - q) * x)^(1 / (1 - q))
@. q_gaussian(x, p) = p[1] * q_exponential(-p[2] * x^2, p[3])

# Example data
xdata = range(0, stop=10, length=20)
ydata = q_gaussian(xdata, [1.0, 0.1, 1.5]) + 0.01 * randn(length(xdata))
p0 = [1.0, 0.1, 1.5]  # Initial guess for parameters

# Fit the q-Gaussian model to the data
fit = curve_fit(q_gaussian, xdata, ydata, p0)

# Optional: Bounds and fit with bounds
#lb = [0.5, 0.0, 1.0]  # Lower bounds
#ub = [1.5, 0.2, 2.0]  # Upper bounds
fit_bounds = curve_fit(q_gaussian, xdata, ydata, p0)#, lower=lb, upper=ub)

# Get standard errors and confidence intervals
sigma = stderror(fit)
margin_of_error = margin_error(fit, 0.05)

# Jacobian function for the q-Gaussian model
function jacobian_model(x, p)
    A, B, q = p
    J = Array{Float64}(undef, length(x), length(p))
    @. J[:, 1] = q_exponential(-B * x^2, q)
    @. J[:, 2] = -2 * B * x^2 * A * q * q_exponential(-B * x^2, q)
    @. J[:, 3] = A * q_exponential(-B * x^2, q) * log(1 + (1 - q) * (-B * x^2))
    J
end

fit= curve_fit(q_gaussian, jacobian_model, xdata, ydata, p0)

