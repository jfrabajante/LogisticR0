# Set working directory
setwd("/Users/jomarrabajante/Downloads/sir/simulate")  # Don't forget to change this

# Install necessary packages (run this section only once)
# install.packages("tidyverse")
# install.packages("caret")
# install.packages("pROC")
# install.packages("broom")
# install.packages("lhs")
# install.packages("deSolve")
# install.packages("sandwich")
# install.packages("lmtest")
# install.packages("plotly")

# Load necessary libraries
library(tidyverse)
library(caret)
library(pROC)
library(broom)
library(lhs)
library(deSolve)
library(sandwich)
library(lmtest)
library(plotly)

##################### EPIDEMIC MODEL SIMULATIONS

# Compartmental model function using deSolve 
seir_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    #N <- S + I + R  # Frequency dependent
    N <- 1  # Density dependent
    
    incidence <- min(beta * S * I / N, S)
    deltaS <- -incidence
    deltaI <- incidence - gamma * I
    deltaR <- gamma * I
    
    # Assure non-negativity
    dS <- max(deltaS, -S)
    dI <- max(deltaI, -I)
    dR <- max(deltaR, -R)
    
    dZ <- ifelse (time > 1/gamma, 0, ifelse (I > 0, incidence/I, 0))  # To record average secondary cases infected by the index case (approximate basic reproduction number R_0) - be careful in defining this as this is model-based
    
    dC <- incidence  # Cumulative cases (to check approximate final epidemic size)
    
    list(c(dS, dI, dR, dZ, dC))
  })
}

# Perform Latin Hypercube Sampling
set.seed(42)

# Define the parameter ranges
param_ranges <- data.frame(
  beta = c(0, 1),
  gamma = c(0, 1)
)

# Generate samples using Latin Hypercube Sampling
n_samples <- 2*400
samples <- randomLHS(n_samples, ncol(param_ranges))
colnames(samples) <- colnames(param_ranges)
samples <- as.data.frame(samples)

# Scale samples to parameter ranges
for (i in 1:ncol(samples)) {
  param_min <- param_ranges[1, i]  # Minimum value for parameter i
  param_max <- param_ranges[2, i]  # Maximum value for parameter i
  samples[, i] <- samples[, i] * (param_max - param_min) + param_min
}

# Define initial state and time steps
initial_state <- c(S = 99/100, I = 1/100, R = 0, Z = 0, C = 1/100)  
time_steps <- seq(0, 500, by = 1)

# Initialize a dataframe to store results
results <- data.frame()

# Run the Simulation Multiple Times
for (i in 1:n_samples) {
  params <- as.list(samples[i, ])
  sim_results <- ode(y = initial_state, times = time_steps, func = seir_model, parms = params)
  sim_results <- as.data.frame(sim_results)
  
  # Record the maximum value of I, Z and C
  max_I <- max(sim_results$I)
  max_Z <- max(sim_results$Z)
  max_C <- max(sim_results$C)
  
  # Check for outbreak (I > threshold)
  outbreak_thresholdI <- 1/100  # Define an appropriate threshold for an outbreak 
  outbreak_indicatorI <- ifelse(max_I > outbreak_thresholdI, 1, 0)
  
  outbreak_threshold <- 1  # Define an appropriate threshold for an outbreak R0>1
  outbreak_indicator <- ifelse(max_Z > outbreak_threshold, 1, 0)
  
  outbreak_threshold2 <- 2  # Define an appropriate threshold for an outbreak R0>2
  outbreak_indicator2 <- ifelse(max_Z > outbreak_threshold2, 1, 0)
  
  outbreak_threshold3 <- 3  # Define an appropriate threshold for an outbreak R0>3
  outbreak_indicator3 <- ifelse(max_Z > outbreak_threshold3, 1, 0)
  
  outbreak_threshold4 <- 4  # Define an appropriate threshold for an outbreak R0>4 
  outbreak_indicator4 <- ifelse(max_Z > outbreak_threshold4, 1, 0)
  
  outbreak_threshold5 <- 5  # Define an appropriate threshold for an outbreak R0>5
  outbreak_indicator5 <- ifelse(max_Z > outbreak_threshold5, 1, 0)
  
  outbreak_threshold6 <- 6  # Define an appropriate threshold for an outbreak R0>6
  outbreak_indicator6 <- ifelse(max_Z > outbreak_threshold6, 1, 0)
  
  outbreak_threshold7 <- 7  # Define an appropriate threshold for an outbreak R0>7
  outbreak_indicator7 <- ifelse(max_Z > outbreak_threshold7, 1, 0)
    
  # Store results
  results <- rbind(results, cbind(samples[i, ], max_Z = max_Z, outbreak = outbreak_indicator, outbreak2 = outbreak_indicator2, outbreak3 = outbreak_indicator3, outbreak4 = outbreak_indicator4, outbreak5 = outbreak_indicator5, outbreak6 = outbreak_indicator6, outbreak7 = outbreak_indicator7, max_I = max_I, outbreakI = outbreak_indicatorI, total_cases = max_C))
}

# Save Parameter Sets, Outbreak Indicators, and Max Values
write.csv(results, file = "simulation_results.csv", row.names = FALSE)

##################### MACHINE LEARNING PART

# Load the simulation results
data <- read.csv("simulation_results.csv")

# Ensure all values are numeric and there are no NAs
data <- data %>%
  mutate(across(everything(), as.numeric)) %>%
  drop_na()

##################### Outbreak R0>1

# Fit the logistic regression model
model <- glm(outbreak ~ beta + gamma, data = data, family = binomial)  # Change outbreak to fit other indicators

# Calculate robust standard errors
robust_se <- sqrt(diag(vcovHC(model, type = "HC0")))

# Extract coefficients and p-values with robust standard errors
coefs <- summary(model)$coefficients
robust_p_values <- coeftest(model, vcov. = vcovHC(model, type = "HC0"))[, "Pr(>|z|)"]

# Save coefficients with robust standard errors
robust_results <- data.frame(
  Term = rownames(coefs),
  Estimate = coefs[, "Estimate"],
  Std_Error = robust_se,
  z_value = coefs[, "Estimate"] / robust_se,
  Pr_z = robust_p_values
)
write.csv(robust_results, file = "robust_model_results.csv", row.names = FALSE)

# Calculate the predicted probabilities
data$predicted_prob <- predict(model, type = "response")

# Calculate the predicted classes based on a threshold of 0.5
data$predicted_class <- ifelse(data$predicted_prob > 0.5, 1, 0)

# Confusion matrix
conf_matrix <- table(data$outbreak, data$predicted_class)  # Change outbreak to fit other indicators
write.csv(as.data.frame(conf_matrix), file = "confusion_matrix.csv", row.names = FALSE)

# Performance metrics
confusion <- confusionMatrix(as.factor(data$predicted_class), as.factor(data$outbreak))  # Change outbreak to fit other indicators
precision <- confusion$byClass["Pos Pred Value"]
sensitivity <- confusion$byClass["Sensitivity"]
specificity <- confusion$byClass["Specificity"]
auc_value <- auc(roc(data$outbreak, data$predicted_prob))  # Change outbreak to fit other indicators
mse <- mean((data$predicted_prob - data$outbreak)^2)  # Change outbreak to fit other indicators

performance_metrics <- data.frame(
  Metric = c("Accuracy", "AUC", "Precision", "Sensitivity", "Specificity", "H0 vs H1 p-value", "MSE"),
  Value = c(
    mean(data$predicted_class == data$outbreak),  # Change outbreak to fit other indicators
    auc_value,
    precision,
    sensitivity,
    specificity,
    min(robust_p_values),
    mse
  )
)
write.csv(performance_metrics, file = "performance_metrics.csv", row.names = FALSE)

print(paste("Accuracy:", mean(data$predicted_class == data$outbreak)))  # Change outbreak to fit other indicators
print(paste("Precision:", precision))
print(paste("Sensitivity:", sensitivity))
print(paste("Specificity:", specificity))
print(paste("H0 vs H1 p-value:", min(robust_p_values)))
print(paste("Mean Squared Error:", mse))

# Plot ROC Curve
roc_curve <- roc(data$outbreak, data$predicted_prob)  # Change outbreak to fit other indicators
plot(roc_curve, main = "ROC Curve", col = "blue")
auc(roc_curve)  # Print AUC value

# Define the logistic regression equation based on the model coefficients
intercept <- coef(model)["(Intercept)"]
beta_coeff <- coef(model)["beta"]
gamma_coeff <- coef(model)["gamma"]

# Print the logistic regression equation
equation <- paste0("logit(p) = ", round(intercept, 3), 
                   " + ", round(beta_coeff, 3), "*beta + ", 
                   round(gamma_coeff, 3), "*gamma")

# Display the logistic regression equation
print(paste("Logistic Regression Equation:"))
print(equation)

# Convert the logit function to probability
equation_probability <- paste0("p = 1 / (1 + exp(-(", round(intercept, 3), 
                               " + ", round(beta_coeff, 3), "*beta + ", 
                               round(gamma_coeff, 3), "*gamma)))")
print(paste("Logistic Regression Equation (Probability Form):"))
print(equation_probability)  # This is the R_0-like epidemic threshold

data$predicted_prob <- 1 / (1 + exp(-(intercept + beta_coeff * data$beta + gamma_coeff * data$gamma)))

outdiff <- abs(data$outbreak - data$outbreakI)
sum(outdiff)

# Save updated data with predicted probabilities
write.csv(data, file = "updated_simulation_results.csv", row.names = FALSE)

# Plot 3D surface of beta, gamma, and predicted probability
fig <- plot_ly(data, x = ~beta, y = ~gamma, z = ~predicted_prob, type = "scatter3d", mode = "markers", 
               marker = list(size = 3, color = ~predicted_prob, colorscale = "Viridis", showscale = TRUE))

# Add a surface plot to show how `predicted_prob` changes with `beta` and `gamma`
fig <- fig %>% add_surface(
  x = seq(min(data$beta), max(data$beta), length.out = 50),
  y = seq(min(data$gamma), max(data$gamma), length.out = 50),
  z = outer(seq(min(data$beta), max(data$beta), length.out = 50), 
            seq(min(data$gamma), max(data$gamma), length.out = 50), 
            function(b, g) 1 / (1 + exp(-(intercept + beta_coeff * b + gamma_coeff * g))))
)

fig <- fig %>% layout(scene = list(xaxis = list(title = 'Beta'),
                                   yaxis = list(title = 'Gamma'),
                                   zaxis = list(title = 'Probability')),
                      title = '3D Plot of Beta, Gamma, and Predicted Probability')

fig

##################### Outbreak R0>2

# Fit the logistic regression model
model <- glm(outbreak2 ~ beta + gamma, data = data, family = binomial)  # Change outbreak to fit other indicators

# Calculate robust standard errors
robust_se <- sqrt(diag(vcovHC(model, type = "HC0")))

# Extract coefficients and p-values with robust standard errors
coefs <- summary(model)$coefficients
robust_p_values <- coeftest(model, vcov. = vcovHC(model, type = "HC0"))[, "Pr(>|z|)"]

# Save coefficients with robust standard errors
robust_results <- data.frame(
  Term = rownames(coefs),
  Estimate = coefs[, "Estimate"],
  Std_Error = robust_se,
  z_value = coefs[, "Estimate"] / robust_se,
  Pr_z = robust_p_values
)
write.csv(robust_results, file = "robust_model_results2.csv", row.names = FALSE)

# Calculate the predicted probabilities
data$predicted_prob <- predict(model, type = "response")

# Calculate the predicted classes based on a threshold of 0.5
data$predicted_class <- ifelse(data$predicted_prob > 0.5, 1, 0)

# Confusion matrix
conf_matrix <- table(data$outbreak2, data$predicted_class)  # Change outbreak to fit other indicators
write.csv(as.data.frame(conf_matrix), file = "confusion_matrix2.csv", row.names = FALSE)

# Performance metrics
confusion <- confusionMatrix(as.factor(data$predicted_class), as.factor(data$outbreak2))  # Change outbreak to fit other indicators
precision <- confusion$byClass["Pos Pred Value"]
sensitivity <- confusion$byClass["Sensitivity"]
specificity <- confusion$byClass["Specificity"]
auc_value <- auc(roc(data$outbreak2, data$predicted_prob))  # Change outbreak to fit other indicators
mse <- mean((data$predicted_prob - data$outbreak2)^2)  # Change outbreak to fit other indicators

performance_metrics <- data.frame(
  Metric = c("Accuracy", "AUC", "Precision", "Sensitivity", "Specificity", "H0 vs H1 p-value", "MSE"),
  Value = c(
    mean(data$predicted_class == data$outbreak2),  # Change outbreak to fit other indicators
    auc_value,
    precision,
    sensitivity,
    specificity,
    min(robust_p_values),
    mse
  )
)
write.csv(performance_metrics, file = "performance_metrics2.csv", row.names = FALSE)

print(paste("Accuracy:", mean(data$predicted_class == data$outbreak2)))  # Change outbreak to fit other indicators
print(paste("Precision:", precision))
print(paste("Sensitivity:", sensitivity))
print(paste("Specificity:", specificity))
print(paste("H0 vs H1 p-value:", min(robust_p_values)))
print(paste("Mean Squared Error:", mse))

# Plot ROC Curve
roc_curve <- roc(data$outbreak2, data$predicted_prob)  # Change outbreak to fit other indicators
plot(roc_curve, main = "ROC Curve", col = "blue")
auc(roc_curve)  # Print AUC value

# Define the logistic regression equation based on the model coefficients
intercept <- coef(model)["(Intercept)"]
beta_coeff <- coef(model)["beta"]
gamma_coeff <- coef(model)["gamma"]

# Print the logistic regression equation
equation <- paste0("logit(p) = ", round(intercept, 3), 
                   " + ", round(beta_coeff, 3), "*beta + ", 
                   round(gamma_coeff, 3), "*gamma")

# Display the logistic regression equation
print(paste("Logistic Regression Equation (R0>2):"))
print(equation)

# Convert the logit function to probability
equation_probability <- paste0("p = 1 / (1 + exp(-(", round(intercept, 3), 
                               " + ", round(beta_coeff, 3), "*beta + ", 
                               round(gamma_coeff, 3), "*gamma)))")
print(paste("Logistic Regression Equation (Probability Form; R0>2):"))
print(equation_probability)  # This is the R_0-like epidemic threshold

data$predicted_prob <- 1 / (1 + exp(-(intercept + beta_coeff * data$beta + gamma_coeff * data$gamma)))

# Save updated data with predicted probabilities
write.csv(data, file = "updated_simulation_results2.csv", row.names = FALSE)

# Plot 3D surface of beta, gamma, and predicted probability
fig <- plot_ly(data, x = ~beta, y = ~gamma, z = ~predicted_prob, type = "scatter3d", mode = "markers", 
               marker = list(size = 3, color = ~predicted_prob, colorscale = "Viridis", showscale = TRUE))

# Add a surface plot to show how `predicted_prob` changes with `beta` and `gamma`
fig <- fig %>% add_surface(
  x = seq(min(data$beta), max(data$beta), length.out = 50),
  y = seq(min(data$gamma), max(data$gamma), length.out = 50),
  z = outer(seq(min(data$beta), max(data$beta), length.out = 50), 
            seq(min(data$gamma), max(data$gamma), length.out = 50), 
            function(b, g) 1 / (1 + exp(-(intercept + beta_coeff * b + gamma_coeff * g))))
)

fig <- fig %>% layout(scene = list(xaxis = list(title = 'Beta'),
                                   yaxis = list(title = 'Gamma'),
                                   zaxis = list(title = 'Probability')),
                      title = '3D Plot of Beta, Gamma, and Predicted Probability (R0>2)')

fig

##################### Outbreak R0>3

# Fit the logistic regression model
model <- glm(outbreak3 ~ beta + gamma, data = data, family = binomial)  # Change outbreak to fit other indicators

# Calculate robust standard errors
robust_se <- sqrt(diag(vcovHC(model, type = "HC0")))

# Extract coefficients and p-values with robust standard errors
coefs <- summary(model)$coefficients
robust_p_values <- coeftest(model, vcov. = vcovHC(model, type = "HC0"))[, "Pr(>|z|)"]

# Save coefficients with robust standard errors
robust_results <- data.frame(
  Term = rownames(coefs),
  Estimate = coefs[, "Estimate"],
  Std_Error = robust_se,
  z_value = coefs[, "Estimate"] / robust_se,
  Pr_z = robust_p_values
)
write.csv(robust_results, file = "robust_model_results3.csv", row.names = FALSE)

# Calculate the predicted probabilities
data$predicted_prob <- predict(model, type = "response")

# Calculate the predicted classes based on a threshold of 0.5
data$predicted_class <- ifelse(data$predicted_prob > 0.5, 1, 0)

# Confusion matrix
conf_matrix <- table(data$outbreak3, data$predicted_class)  # Change outbreak to fit other indicators
write.csv(as.data.frame(conf_matrix), file = "confusion_matrix3.csv", row.names = FALSE)

# Performance metrics
confusion <- confusionMatrix(as.factor(data$predicted_class), as.factor(data$outbreak3))  # Change outbreak to fit other indicators
precision <- confusion$byClass["Pos Pred Value"]
sensitivity <- confusion$byClass["Sensitivity"]
specificity <- confusion$byClass["Specificity"]
auc_value <- auc(roc(data$outbreak3, data$predicted_prob))  # Change outbreak to fit other indicators
mse <- mean((data$predicted_prob - data$outbreak3)^2)  # Change outbreak to fit other indicators

performance_metrics <- data.frame(
  Metric = c("Accuracy", "AUC", "Precision", "Sensitivity", "Specificity", "H0 vs H1 p-value", "MSE"),
  Value = c(
    mean(data$predicted_class == data$outbreak3),  # Change outbreak to fit other indicators
    auc_value,
    precision,
    sensitivity,
    specificity,
    min(robust_p_values),
    mse
  )
)
write.csv(performance_metrics, file = "performance_metrics3.csv", row.names = FALSE)

print(paste("Accuracy:", mean(data$predicted_class == data$outbreak3)))  # Change outbreak to fit other indicators
print(paste("Precision:", precision))
print(paste("Sensitivity:", sensitivity))
print(paste("Specificity:", specificity))
print(paste("H0 vs H1 p-value:", min(robust_p_values)))
print(paste("Mean Squared Error:", mse))

# Plot ROC Curve
roc_curve <- roc(data$outbreak3, data$predicted_prob)  # Change outbreak to fit other indicators
plot(roc_curve, main = "ROC Curve", col = "blue")
auc(roc_curve)  # Print AUC value

# Define the logistic regression equation based on the model coefficients
intercept <- coef(model)["(Intercept)"]
beta_coeff <- coef(model)["beta"]
gamma_coeff <- coef(model)["gamma"]

# Print the logistic regression equation
equation <- paste0("logit(p) = ", round(intercept, 3), 
                   " + ", round(beta_coeff, 3), "*beta + ", 
                   round(gamma_coeff, 3), "*gamma")

# Display the logistic regression equation
print(paste("Logistic Regression Equation (R0>3):"))
print(equation)

# Convert the logit function to probability
equation_probability <- paste0("p = 1 / (1 + exp(-(", round(intercept, 3), 
                               " + ", round(beta_coeff, 3), "*beta + ", 
                               round(gamma_coeff, 3), "*gamma)))")
print(paste("Logistic Regression Equation (Probability Form; R0>3):"))
print(equation_probability)  # This is the R_0-like epidemic threshold

data$predicted_prob <- 1 / (1 + exp(-(intercept + beta_coeff * data$beta + gamma_coeff * data$gamma)))

# Save updated data with predicted probabilities
write.csv(data, file = "updated_simulation_results3.csv", row.names = FALSE)

# Plot 3D surface of beta, gamma, and predicted probability
fig <- plot_ly(data, x = ~beta, y = ~gamma, z = ~predicted_prob, type = "scatter3d", mode = "markers", 
               marker = list(size = 3, color = ~predicted_prob, colorscale = "Viridis", showscale = TRUE))

# Add a surface plot to show how `predicted_prob` changes with `beta` and `gamma`
fig <- fig %>% add_surface(
  x = seq(min(data$beta), max(data$beta), length.out = 50),
  y = seq(min(data$gamma), max(data$gamma), length.out = 50),
  z = outer(seq(min(data$beta), max(data$beta), length.out = 50), 
            seq(min(data$gamma), max(data$gamma), length.out = 50), 
            function(b, g) 1 / (1 + exp(-(intercept + beta_coeff * b + gamma_coeff * g))))
)

fig <- fig %>% layout(scene = list(xaxis = list(title = 'Beta'),
                                   yaxis = list(title = 'Gamma'),
                                   zaxis = list(title = 'Probability')),
                      title = '3D Plot of Beta, Gamma, and Predicted Probability (R0>3)')

fig