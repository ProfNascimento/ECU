## LOADING LIBRARY
library(tidyverse)
library(forecast)
library(bnlearn)
library(qgraph)
library(Metrics)
library(rbmn)
library(cowplot)

set.seed(20)

#######################################################
# DATA IMPORT
novoscasos = read.csv("https://raw.githubusercontent.com/ProfNascimento/spatialBN/main/COVID.csv")

# DATA WRANGLING
DB = novoscasos %>% 
  filter(date <= '2022-02-15') %>% mutate(diag = as.Date(novoscasos$date),
                                                            weekno = as.numeric(diag - first(diag)),
                                                            weekno = (weekno %/% 7) +1) %>% 
  mutate(area = ifelse(state == 'AC' | state == "AM" | state ==  "RR" | state ==  "RO" | state ==  "AP" | state ==  "PA" | state ==  "TO", 'North', 
                       ifelse(state == "MA" | state ==  "PI" | state ==  "CE" | state ==  "RN" | state ==  "PB" | state ==  "PE" | state ==  "AL" | state ==  "SE" | state ==  "BA", 'Northeast',
                              ifelse(state == "PR" | state == "SC" | state ==  "RS", 'South',
                                     ifelse(state == "MG" | state == "SP" | state == "ES" | state == "RJ", 'Southeast', 'Midwest')))) ) %>% 
  group_by(weekno, area) %>%  
  summarise(value = sum(newCases)/10000)

# VISUALIZATION PER REGION
ggplot(DB, aes(x = weekno, y = value, fill = area, colour = area)) +
  geom_line(size = 0.8) +
  scale_fill_manual(values = c(`Midwest` = "#FF1100", Northeast = "#F4E901", North = "#075500", 
                               Southeast = "#FF4A90", South = "#2C13EB")) +
  scale_color_manual(values = c(`Midwest` = "#FF1100", Northeast = "#F4E901", North = "#075500",
                                Southeast = "#FF4A90", South = "#2C13EB")) +
  labs(x = "Epidemiological Week", 
       y = "New cases (per 10k)", fill = "Region:", color = "Region:") +
  theme_minimal() +
  theme(legend.position = "bottom",axis.title.y = element_text(size = 18L), 
        axis.title.x = element_text(size = 18L), axis.text = element_text(size = rel(1.5)),
        legend.text = element_text(size = 16L), legend.title = element_text(size = 16L))

######################################################################
## STEP 01 - SPATIAL (NETWORK) ESTIMATION
DB <- tidyr::pivot_wider(DB, names_from = "area", values_from = "value")
DB[is.na(DB)] <- 0
train <- DB[1:82,-1]
teste <- DB[-c(1:82),-1]

# TABU MODEL (NON-INFORMATIVE SPATIAL PRIOR)
modtabu_train <- boot.strength(train, R = 1000, algorithm = 'tabu', debug = FALSE)
modtabu_train2 = modtabu_train

qgraph(modtabu_train, asize=5,
       legend.cex=0.5,
       edge.color="black",
       color="#80deea", 
       edge.labels=T)

# HEATMAP BN PLOT
modtabu_train %>% mutate(from = fct_relevel(from, "Midwest", "Northeast", "North", "Southeast", "South"),
                         to = fct_relevel(to, "Midwest", "Northeast", "North", "Southeast", "South") ) %>% 
  ggplot( aes(from, to) ) +
  geom_tile( aes(fill = strength), color = "white" ) + scale_fill_gradient(low = "white", high = "blue") +
  ylab(" ") +  xlab(" ") + labs(fill = " ") +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        plot.title = element_text(size=18),
        axis.title = element_text(size=16,face="bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text = element_text(size = rel(1.5)))

modtabu_train_medio <- averaged.network(modtabu_train, threshold = 0.5)
modelstring(modtabu_train_medio)

graphviz.plot(bnlearn::model2network(modelstring(modtabu_train_medio)), layout = "fdp")

# FIT MODEL TABU (NON-INFORMATIVE PRIOR)
modelstring(modtabu_train_medio) 
mod.Tabu = bn.fit(modtabu_train_medio, data = train)

######################################################################
## STEP 02 - TIME SERIES MODELING
# Define the function
results = list()
analyze_model <- function(DB, model) {
  # Step 1: Iterate over each node in the model
  for (list_name in names(model)) {
    # Calculate residuals for the current node
    ress <- DB[[list_name]] - predict(model, node = list_name, data = DB)
    
    # Step 2: Fit ARIMA model to the residuals
    fit <- auto.arima(ress)
    
    # Store the results in the list
    results[[list_name]] <- list(fit = fit)
  }
  # Return the list of results
  return(results)
}

## model = {mod.Tabu, mod.Rod, mod.Rod2, mod.Air, mod.Air2}
( BN.ARIMA = analyze_model(DB, model=mod.Rod) )

# UNIVARIATE RESIDUAL ANALYSIS (BN.ARIMA)
cpgram(residuals(BN.ARIMA$North$fit))
acf(residuals(BN.ARIMA$North$fit))
Box.test(residuals(BN.ARIMA$North$fit), lag = 20, type = "Ljung-Box")
shapiro.test(residuals(BN.ARIMA$North$fit))

# VISUAL COMPARISON MODELS (OBS SERIES, BN FORECAST, BN+ARIMA FORECAST)
FIT_plot <- function(DB, BN.model, list_name) {
  # Step 1: Extract actual values
  actual_values <- DB[[list_name]]
  
  # Step 2: Get predicted values from the Bayesian Network model
  bn_predictions <- predict(BN.model, node = list_name, data = DB)
  
  # Step 3: Get fitted values from the ARIMA model
  ress <- DB[[list_name]] - c(predict(BN.model, node = list_name, data = DB))
  fit <- auto.arima(ress)
  arima_fitted <- fitted(fit)
  
  # Step 4: Extract residuals from the ARIMA model
  residuals_node <- residuals(BN.ARIMA[[list_name]]$fit)
  
  # Create a data frame for plotting
  plot_data <- data.frame(
    Time = 1:length(actual_values),
    Actual = actual_values,
    BN_Predicted = bn_predictions,
    ARIMA_Fitted = arima_fitted + bn_predictions
  )
  
  # Step 5: Plot the actual vs predicted values
  time_series_plot <- ggplot(plot_data, aes(x = Time)) +
    geom_line(aes(y = Actual, color = "Actual"), size = 1) +  # Solid line for actual values
    geom_line(aes(y = BN_Predicted, color = "BN Predicted"), linetype = "dashed", size = 1) +  # Dashed blue line for BN predictions
    geom_line(aes(y = ARIMA_Fitted, color = "BN-ARIMA Fitted"), linetype = "dashed", size = 1) +  # Dashed orange line for ARIMA fitted values
    scale_color_manual(values = c("Actual" = "black", "BN Predicted" = "blue", "BN-ARIMA Fitted" = "orange")) +  # Correct color mapping
    labs(title = paste("Actual vs Predicted -", list_name),
         x = "Time",
         y = "Values",
         color = "Legend") +  # Label for the legend
    theme_minimal() +
    theme(legend.position = c(0.1, 0.9),  # Set legend position to top-left
          legend.background = element_blank(),  # No background for the legend
          legend.key = element_blank())  # Remove any boxes around the legend keys
  
  # Step 6: Create cpgram (Cumulative Periodogram) plot for residuals
  cpgram_plot <- ggplot(data = data.frame(residuals = residuals_node), aes(x = residuals)) +
    stat_ecdf() +
    labs(title = paste(list_name),
         x = "Residuals",
         y = "Cumulative Probability") +
    theme_minimal()
  
  # Step 7: Create acf (Autocorrelation Function) plot for residuals
  acf_values <- acf(residuals_node, plot = FALSE)
  acf_data <- data.frame(lag = acf_values$lag, acf = acf_values$acf)
  
  acf_plot <- ggplot(acf_data, aes(x = lag, y = acf)) +
    geom_bar(stat = "identity") +
    labs(title = paste(list_name),
         x = "Lag",
         y = "ACF") +
    theme_minimal()
  
  # Combine the plots
  print( plot_grid(time_series_plot,
                   plot_grid(cpgram_plot, acf_plot, ncol = 1),ncol=2, rel_widths = c(2,1), scale = c(1,.9)) )
  
  # Step 8: Calculate RMSE for BN-ARIMA fitted values and BN scores
  bn.arima_rmse <- round(rmse(plot_data$Actual, plot_data$ARIMA_Fitted),2)
  bn_rmse <- round(rmse(plot_data$Actual, plot_data$BN_Predicted),2)
  bn_aic <- round(score(bn.net(BN.model), teste, type = "aic-g"),2)
  bn_bic <- round(score(bn.net(BN.model), teste, type = "bic-g"),2)
  bn_bge <- round(score(bn.net(BN.model), teste, type = "bge"),2)
  
  # Step 10: Print the RMSE and BN scores
  print(paste("RMSE (BN-ARIMA):", bn.arima_rmse))
  print(paste("RMSE (BN):", bn_rmse))
  print(paste("BN AIC:", bn_aic))
  print(paste("BN BIC:", bn_bic))
  print(paste("BN BGE:", bn_bge))
}

## list_name = {"North", "Northeast", "Midwest", "Southeast", "South"}
FIT_plot(DB, BN.model = mod.Tabu, 
         list_name = "South")

## REGRESSION COEF HEATMAP --ROAD DIST. MODEL--
bnfit2nbn(mod.Tabu) %>% 
  rmatrix4nbn(stdev = FALSE) %>% 
  reshape2::melt(id.vars = "to") %>%
  mutate(from = fct_relevel(from, "Midwest", "Northeast", "North", "Southeast", "South"),
         to = fct_relevel(to, "Midwest", "Northeast", "North", "Southeast", "South") ) %>% 
  ggplot( aes(from, to) ) +
  geom_tile(aes(fill = round(value,2) )) +
  scale_fill_gradient2(low="darkred", high="darkgreen", guide="colorbar")+
  guides(fill = guide_colourbar(title = "REG. COEF."))+
  ylab("Child Node") + xlab("Parent Node") +
  theme(legend.position = "right",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 18),
        plot.title = element_text(size=18),
        axis.title=element_text(size=18,face="bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text = element_text(size = rel(1.5))) +
  labs(fill = " ") 
