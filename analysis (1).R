

# Q2: ODE model predicts cisplatin sensitivity in NSCLC cell lines
# Q5: ODE model vs ML approaches



#Libraries
library(deSolve)
library(ggplot2)
library(glmnet)
library(randomForest)
library(neuralnet)
library(tidyverse)
library(readxl)
library(gridExtra)

set.seed(42)

#p53 - Mdm2 ODE Model
p53_model <- function(t, y, params) {
  p53  <- y[1]
  mRNA <- y[2]
  Mdm2 <- y[3]

  with(as.list(params), {
    dp53  <- alpha_p - beta_p * p53 - gamma_p * Mdm2 * p53 / (1 + kappa * DS)
    dmRNA <- alpha_m * p53 / (K + p53) - beta_m * mRNA
    dMdm2 <- alpha_M * mRNA - (beta_M + kappa_d * DS) * Mdm2
    return(list(c(dp53, dmRNA, dMdm2)))
  })
}

# Default parameters
default_params <- c(
  alpha_p = 0.9, beta_p  = 0.1, gamma_p = 1.8,
  kappa   = 0.9, alpha_m = 1.5, beta_m  = 0.8,
  K       = 0.5, alpha_M = 0.9, beta_M  = 0.4,
  kappa_d = 1.2, DS      = 0.8
)

# Initial conditions and time
y0    <- c(p53 = 0.5, mRNA = 0.5, Mdm2 = 0.5)
times <- seq(0, 72, by = 0.1)

# Map cell line genomic features to ODE parameters
get_params <- function(row) {
  p <- default_params

  if (row$tp53_mut == 1) {
    p["alpha_m"] <- p["alpha_m"] * 0.25
    p["gamma_p"] <- p["gamma_p"] * 0.50
    p["kappa"]   <- p["kappa"]   * 0.40
  }

  p["alpha_M"] <- p["alpha_M"] * 2^(row$mdm2_expr - 6.0) / 4

  if (row$atm_mut == 1) {
    p["kappa_d"] <- p["kappa_d"] * 0.15
    p["kappa"]   <- p["kappa"]   * 0.20
  }

  return(p)
}

# Extract features from p53 
get_features <- function(p53_trace, times) {
  peak_p53        <- max(p53_trace)
  mean_p53        <- mean(p53_trace)
  apoptosis_score <- mean(pmax(0, p53_trace - 1.2))
  late_p53        <- mean(tail(p53_trace, 240))
  return(c(peak_p53, mean_p53, apoptosis_score, late_p53))
}

# GDSC2 drug sensitivity
gdsc <- read_excel("gdsc2.xlsx")
lung_gdsc <- gdsc[gdsc$TCGA_DESC %in% c("LUAD", "LUSC") &
                    gdsc$DRUG_NAME == "Cisplatin", ]

# DepMap expression
expr <- read.csv("OmicsExpressionProteinCodingGenesTPMLogp1.csv",
                 row.names = 1, check.names = FALSE)

expr_small <- data.frame(
  ModelID   = rownames(expr),
  mdm2_expr = expr[, "MDM2 (4193)"],
  atm_expr  = expr[, "ATM (472)"]
)

# DepMap mutations
muts <- read.csv("OmicsSomaticMutationsMatrixDamaging.csv",
                 row.names = 1, check.names = FALSE)

muts_small <- data.frame(
  ModelID  = muts$ModelID,
  tp53_mut = ifelse(muts[, "TP53 (7157)"] > 0, 1, 0),
  atm_mut  = ifelse(muts[, "ATM (472)"]   > 0, 1, 0)
)

# Model info
model_info <- read.csv("Model.csv")

# Add cell line name
expr_small <- merge(expr_small,
                    model_info[, c("ModelID", "CellLineName")],
                    by = "ModelID")

muts_small <- merge(muts_small,
                    model_info[, c("ModelID", "CellLineName")],
                    by = "ModelID")

# Merge
real_data <- merge(lung_gdsc[, c("CELL_LINE_NAME", "LN_IC50")],
                   expr_small[, c("CellLineName", "mdm2_expr", "atm_expr")],
                   by.x = "CELL_LINE_NAME", by.y = "CellLineName")

real_data <- merge(real_data,
                   muts_small[, c("CellLineName", "tp53_mut", "atm_mut")],
                   by.x = "CELL_LINE_NAME", by.y = "CellLineName")

real_data <- unique(real_data)

cat("Cell lines after merging:", nrow(real_data), "\n")

#ODE Simulation
n <- nrow(real_data)
ode_feat <- matrix(NA, nrow = n, ncol = 4)
colnames(ode_feat) <- c("peak_p53", "mean_p53", "apoptosis_score", "late_p53")

for (i in 1:n) {
  p         <- get_params(real_data[i, ])
  out       <- ode(y = y0, times = times, func = p53_model, parms = p)
  ode_feat[i, ] <- get_features(as.data.frame(out)$p53, times)
}

ode_feat <- as.data.frame(ode_feat)

#Q2: ODE Features → Predict IC50 
X_ode        <- as.matrix(ode_feat)
y            <- real_data$LN_IC50
X_ode_scaled <- scale(X_ode)

cv_ridge <- cv.glmnet(X_ode_scaled, y, alpha = 0)
pred_ode <- as.vector(predict(cv_ridge, X_ode_scaled, s = "lambda.min"))
cor_ode  <- cor(y, pred_ode)

cat("Q2 - ODE + Ridge:  Pearson r =", round(cor_ode, 3), "\n")

#Q5: ML Models on Raw Genomic Features
X_gen        <- as.matrix(real_data[, c("tp53_mut", "mdm2_expr",
                                          "atm_expr", "atm_mut")])
X_gen_scaled <- scale(X_gen)

# Elastic Net
cv_enet   <- cv.glmnet(X_gen_scaled, y, alpha = 0.5)
pred_enet <- as.vector(predict(cv_enet, X_gen_scaled, s = "lambda.min"))
r_enet    <- cor(y, pred_enet)

# Random Forest
rf_model <- randomForest(x = X_gen, y = y, ntree = 200)
pred_rf  <- predict(rf_model, X_gen)
r_rf     <- cor(y, pred_rf)

# Neural Network
nn_data  <- as.data.frame(cbind(X_gen_scaled, LN_IC50 = y))
nn_model <- neuralnet(LN_IC50 ~ tp53_mut + mdm2_expr + atm_expr + atm_mut,
                      data = nn_data, hidden = c(4, 2), linear.output = TRUE)
pred_nn  <- as.vector(predict(nn_model, nn_data))
r_nn     <- cor(y, pred_nn)

# Results table
results <- data.frame(
  Model     = c("ODE + Ridge", "Elastic Net", "Random Forest", "Neural Network"),
  Pearson_r = round(c(cor_ode, r_enet, r_rf, r_nn), 3),
  R_squared = round(c(cor_ode^2, r_enet^2, r_rf^2, r_nn^2), 3)
)

print(results)


#Figures
df_p1 <- data.frame(
  actual    = y,
  predicted = pred_ode,
  tp53      = factor(real_data$tp53_mut, labels = c("WT", "Mutant"))
)

p1 <- ggplot(df_p1, aes(x = actual, y = predicted, colour = tp53)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(title = "Q2: ODE model predictions vs actual cisplatin IC50",
       subtitle = paste("Pearson r =", round(cor_ode, 3), "| n = 57 NSCLC cell lines"),
       x = "Actual ln(IC50)", y = "Predicted ln(IC50)",
       colour = "TP53 status") +
  theme_classic()

p2 <- ggplot(results, aes(x = Model, y = Pearson_r, fill = Model)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  geom_text(aes(label = Pearson_r), vjust = -0.5) +
  ylim(0, 0.85) +
  labs(title = "Q5: ODE vs ML — cisplatin sensitivity prediction",
       subtitle = "Real GDSC2 + DepMap data | NSCLC cell lines",
       x = "", y = "Pearson r") +
  theme_classic()