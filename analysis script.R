# Load necessary libraries
library(metafor)
library(clubSandwich)
library(ggplot2)
library(dplyr)

# Read uploaded data
dat1 <- read.csv('360 sample data.csv')

#calculate ES
data <- escalc(measure="SMD", m1i=Exp_mean, sd1i=Exp_sd, n1i=Exp_n,
               m2i=Ctrl_mean, sd2i=Ctrl_sd, n2i=Ctrl_n, data=dat1)

#display dataset with ES and variance
data

#save .csv file with ES data. This goes into working directory
write.csv(dat1, file = "ESdata.csv")

# Define rho value
rho <- 0.6

# Calculate covariance matrix for rho
V <- vcalc(vi, cluster = Study, obs = ES_number, data = data, rho = rho)

# Compute CHE results for rho
CHEresult <- rma.mv(yi, V,
                    random = ~ 1 | Study/ES_number,
                    method = 'REML',
                    test = 't',
                    dfs = 'contain',
                    data = data,
                    Sparse = TRUE)

# Compute robust CHE results for rho
CHEresultrobust <- robust(CHEresult, cluster = Study, clubSandwich = TRUE, digits = 3)
CHEresultrobust

# Define rho value - 0.2
rho_lower <- 0.4

# Calculate covariance matrix for rho - 0.2
V_lower <- vcalc(vi, cluster = Study, obs = ES_number, data = data, rho = rho_lower)

# Compute CHE results for rho - 0.2
CHEresult_lower <- rma.mv(yi, V_lower,
                          random = ~ 1 | Study/ES_number,
                          method = 'REML',
                          test = 't',
                          dfs = 'contain',
                          data = data,
                          Sparse = TRUE)

# Compute robust CHE results for rho - 0.2
CHEresultrobust_lower <- robust(CHEresult_lower, cluster = Study, clubSandwich = TRUE, digits = 3)
CHEresultrobust_lower

# Define rho value + 0.2
rho_upper <- 0.8

# Calculate covariance matrix for rho + 0.2
V_upper <- vcalc(vi, cluster = Study, obs = ES_number, data = data, rho = rho_upper)

# Compute CHE results for rho + 0.2
CHEresult_upper <- rma.mv(yi, V_upper,
                          random = ~ 1 | Study/ES_number,
                          method = 'REML',
                          test = 't',
                          dfs = 'contain',
                          data = data,
                          Sparse = TRUE)

# Compute robust CHE results for rho + 0.2
CHEresultrobust_upper <- robust(CHEresult_upper, cluster = Study, clubSandwich = TRUE, digits = 3)
CHEresultrobust_upper



#################################################################
##Variance#######################################################

# Teach R i-squared functions
# Please visit the following URL, copy the code, and paste it into your R console and hit enter before proceeding with your code:
# https://raw.githubusercontent.com/MathiasHarrer/dmetar/master/R/mlm.variance.distribution.R

# Calculate I-squared values and variance distribution
i2 <- mlm.variance.distribution(CHEresult)

# Print results and total I2
i2


#################################################################
##Outliers and Influence#########################################

# Compute robust CHE results
robrveinf <- robust(CHEresult, cluster = data$Study, clubSandwich = TRUE, digits = 3)

# Perform outlier analysis
data$upperci <- data$yi + 1.96 * sqrt(data$vi)
data$lowerci <- data$yi - 1.96 * sqrt(data$vi)
data$outlier <- data$upperci < robrveinf$ci.lb | data$lowerci > robrveinf$ci.ub

# Calculate Cook's Distance
cooks <- cooks.distance(robrveinf)

# Calculate DFBETAS
dfbetas <- dfbetas(CHEresult)

#calculate hat values
hatvalues <- hatvalues(robrveinf)

# Calculate p/k
p <- length(coef(robrveinf))
k <- nrow(data)
if (length(coef(robrveinf)) > 1) { p <- p - 1 }
p_over_k <- 3 * (p / k)
hat_flag <- ifelse(hatvalues > p_over_k, 'TRUE', '')

# Combine influence metrics
influence <- data.frame(Study = data$Study,
                        effect_size = data$yi,
                        outlier = ifelse(data$outlier == FALSE, '', 'TRUE'),
                        cooks = cooks,
                        cooks_flag = ifelse(cooks > 0.5, 'TRUE', ''),
                        dfbetas = dfbetas,
                        dfbetas_flag = ifelse(abs(dfbetas) > 1, 'TRUE', ''),
                        hatvalues = hatvalues,
                        hat_flag = hat_flag)

# Adjust column names
colnames(influence)[which(colnames(influence) == 'intrcpt')] <- 'dfbetas*'
colnames(influence)[which(colnames(influence) == 'intrcpt.1')] <- 'dfbetas_flag'

# Add explanatory row
new_row <- data.frame(
  Study = '*Note that DFBETAS is based on the three-level CHE model and not the three-level CHE RVE model as the other metrics are.',
  effect_size = NA,
  outlier = NA,
  cooks = NA,
  cooks_flag = NA,
  dfbetas = NA,
  dfbetas_flag = NA,
  hatvalues = NA,
  hat_flag = NA
)
colnames(new_row) <- colnames(influence)
influence <- rbind(influence, new_row)
influence


#################################################################
##Categorical Moderator Analysis#########################################

# Filter data based on moderator to remove rows where n=1
excluded_levels <- data %>% group_by(Moderator = get('Ctrl_c')) %>% summarise(Count = n(), .groups = 'drop') %>% filter(Count == 1) %>% pull(Moderator)

filtered_data <- data %>% filter(!(get('Ctrl_c') %in% excluded_levels))

# Impute covariance matrix
rho <- 0.60
V_RVE <- vcalc(vi, cluster = Study, obs = ES_number, data = filtered_data, rho = rho)

# Run meta-analysis with intercept for Test of mod
result_tom <- rma.mv(yi = filtered_data$yi, V = V_RVE, mods = ~ factor(filtered_data[['Ctrl_c']]), random = ~ 1 | Study/ES_number, method = 'REML', test = 't', data = filtered_data)
robust_result_tom <- robust(result_tom, cluster = filtered_data$Study, clubSandwich = TRUE, digits = 3)

# Print results for test of moderator
print(robust_result_tom)

# Run meta-analysis  without intercept for table
result <- rma.mv(yi = filtered_data$yi, V = V_RVE, mods = ~ -1 + factor(filtered_data[['Ctrl_c']]), random = ~ 1 | Study/ES_number, method = 'REML', test = 't', data = filtered_data)
robust_result <- robust(result, cluster = filtered_data$Study, clubSandwich = TRUE, digits = 3)

# Print results for table
print(robust_result)



#################################################################
##Plots#########################################

# Run meta-analysis with CHE & robust variance estimation 
result <- rma.mv(yi = data$yi, V = V, random = ~ 1 | Study / ES_number, method = 'REML', test = 't', data = data)
robust_result <- robust(result, cluster = data$Study, clubSandwich = TRUE, digits = 3)

#create comparison level forest plot
forest(robust_result, slab = data$Study, mlab = paste('RE CHE RVE Three-Level Model, rho =', rho), header = TRUE)

# Create comparison-level funnel plot
funnel(robust_result, slab = data$Study, mlab = paste('RE CHE RVE Three-Level Model, rho =', rho), header = TRUE)

# Aggregated data by Study
aggregated_data <- aggregate(cbind(yi, vi) ~ Study, data = data, function(x) mean(x, na.rm = TRUE))

# Impute covariance matrix for aggregated data
V_aggregated <- vcalc(aggregated_data$vi, cluster = aggregated_data$Study, data = aggregated_data, rho = rho)

# Run meta-analysis for aggregated data
result_aggregated <- rma.mv(yi = aggregated_data$yi, V = V_aggregated, random = ~ 1 | Study, method = 'REML', test = 't', data = aggregated_data)
robust_result_aggregated <- robust(result_aggregated, cluster = aggregated_data$Study, clubSandwich = TRUE, digits = 3)

# Generate and save forest plot for aggregated data
forest(robust_result_aggregated, slab = aggregated_data$Study, header = TRUE, annotate = TRUE, addfit = FALSE, mlab = '')