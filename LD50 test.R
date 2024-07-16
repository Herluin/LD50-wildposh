# Herluin Saillard, with Antoine Gekière (Laboratory of Zoology, University of Mons, Belgium) - 2024

#-------------------------   Libraries and data   ----

library(dplyr)
library(drc)
library(ggplot2)
library(ggpubr)
library(patchwork)

data <- file_name


#-------------------------   Dose Response Curve and LC50 dataframe   ----
# OECD (2006), Current approaches in the statistical analysis of ecotoxicity data: a guidance to application
# Ritz et al (2015), Dose-response analysis using R

analyze <- function(data, species_var, stage_var, pesticide_var, exposition_var) {
  
## Antoine Gekière's code DRC - April 2024 ----
  
  # Filter the data for the specific substance and remove NA
  mortality_data <- data %>%
    filter((species == species_var & stage == stage_var & pesticide == pesticide_var | pesticide == "Negative") & exposition == exposition_var) %>% # need different names (eg. "stage = stage" don't work)
    filter(!is.na(censored))
  
  # Generate dataset for analysis
  data_drc <- mortality_data %>%
    group_by(conc_nomil) %>%  # replace "conc_nomil" by the name of the column containing the concentration
    summarise(total = n(), number_dead = sum(censored))
  
  # Calculate the dose-response curve
  drc <- drm(number_dead/total ~ conc_nomil, weights = total, data = data_drc,
             fct = LL.3(fixed = c(NA, 1, NA), names = c("Slope", "Upper Limit", "ED50")), type = "binomial")
  sum_drc <- summary(drc)
  
  # Estimate effective dose ED50
  ED_results <- ED(drc, c(50), interval = "delta")
  
  # LC50 and confidence interval
  LC50 <- ED_results[1]
  standard_error <- ED_results[2]
  confidence_interval <- ED_results[3:4]
  
  # Conversion to µg / bee
  if (exposition_var == "oral") {
    LD50 <- (LC50 / 1000) * 40
    standard_error <- (standard_error / 1000) * 40
    confidence_interval <- (confidence_interval / 1000) * 40
  } else if (exposition_var == "contact") {
    LD50 <- (LC50 / 1000) * 2
    standard_error <- (standard_error / 1000) * 2
    confidence_interval <- (confidence_interval / 1000) * 2
  }

  # Generate new data for the plot
  newdata <- expand.grid(conc_nomil = seq(0, max(data_drc$conc_nomil), by = 0.1))
  pm <- predict(drc, newdata = newdata, interval = "confidence")
  newdata$p <- pm[,1] 
  newdata$pmin <- pm[,2] 
  newdata$pmax <- pm[,3]
  
## Create the dose response curve ----
  plot <- ggplot(data_drc,
                 aes(x = conc_nomil,
                     y = (number_dead/total))) +
    geom_ribbon(data=newdata,
                aes(x=conc_nomil, y=p, ymin=pmin, ymax=pmax),
                alpha=0.2) + 
    geom_line(data=newdata,
              aes(x=conc_nomil, y=p)) +
    geom_point()+
    annotate("segment", x = 0, y = 0.5, xend = LC50, yend = 0.5, color = "blue", linetype = "dashed") +
    annotate("segment", x = LC50, y = 0, xend = LC50, yend = 0.5, color = "blue", linetype = "dashed") +
    annotate("text", x = 0, y = 0.90, label = paste("LC50:", round(LC50, 2), "mg / L"),
             vjust = -1, hjust = 0, size = 4, color = "black", fontface = "bold")+
    annotate("text", x = 0, y = 0.80, label = paste("LD50:", round(LD50, 2), "µg / bee"),
             vjust = -1, hjust = 0, size = 4, color = "black", fontface = "bold")+
    labs(title = paste("DRC for", species_var,  "/", stage_var, "/", exposition_var, "/", pesticide_var),
         x = "Concentration (mg / L)",
         y = "Mortality") +
    theme_bw()
  
## Create dataframe for LC50 and LD50 ----
 # Extract the coefficients and filter for relevant parameters
  coefficients <- as.data.frame(sum_drc$coefficients)
  coefficients <- coefficients %>%
    filter(!grepl("Slope", rownames(coefficients)))
  
  # Create a dataframe with LC50, SE, t-value and p-value
  LC50_data <- data.frame(
    Parameter = paste(species_var, stage_var, pesticide_var, exposition_var),
    LC50_mg_L = coefficients[1, "Estimate"],
    Standard_Error = coefficients[1, "Std. Error"],
    t_value = coefficients[1, "t-value"],
    p_value = coefficients[1, "p-value"])
  
  # Create a dataframe with LD50, SE and confidence interval
  LD50_data <- data.frame(
    Parameter = paste(species_var, stage_var, pesticide_var, exposition_var),
    LD50_µg_per_bee = LD50,
    Standard_Error = standard_error,
    Confidence_interval_lower = confidence_interval[1],
    Confidence_interval_upper = confidence_interval[2])
  
## Var lists ----
  return(list(drc = drc, plot = plot, LD50 = LD50, standard_Error = standard_error, sum_drc = sum_drc, LC50_data = LC50_data, LD50_data = LD50_data))
}

#-------------------------   Modalities   ----

  # use Copilot to generate the code for the other species, pesticides, ...

bombus_adult_acetamiprid_oral <- analyze(data, "bombus", "adult", "Acetamiprid", "oral")
bombus_adult_acetamiprid_contact <- analyze(data, "bombus", "adult", "Acetamiprid", "contact")
bombus_adult_cypermethrin_oral <- analyze(data, "bombus", "adult", "Cypermethrin", "oral")
bombus_adult_cypermethrin_contact <- analyze(data, "bombus", "adult", "Cypermethrin", "contact")

bombus_larvae_acetamiprid_contact <- analyze(data, "bombus", "larvae", "Acetamiprid", "contact")
bombus_larvae_cypermethrin_contact <- analyze(data, "bombus", "larvae", "Cypermethrin", "contact")


#-------------------------   Display drc plot and text table   ----

  # use Copilot to generate the code for the other species, pesticides, ...

# Display the dose response curve
print(bombus_adult_acetamiprid_oral$plot)
print(bombus_adult_acetamiprid_contact$plot)
print(bombus_adult_cypermethrin_oral$plot)
print(bombus_adult_cypermethrin_contact$plot)

print(bombus_larvae_acetamiprid_contact$plot)
print(bombus_larvae_cypermethrin_contact$plot)

# Create and display the text table
LC50_data <- rbind(bombus_adult_acetamiprid_oral$LC50_data, bombus_adult_acetamiprid_contact$LC50_data,
                   bombus_adult_cypermethrin_oral$LC50_data, bombus_adult_cypermethrin_contact$LC50_data,
                   bombus_larvae_acetamiprid_contact$LC50_data,
                   bombus_larvae_cypermethrin_contact$LC50_data)
ggtexttable(LC50_data, rows = NULL, theme = ttheme("lBlueWhite", base_size = 10))

LD50_data <- rbind(bombus_adult_acetamiprid_oral$LD50_data, bombus_adult_acetamiprid_contact$LD50_data,
                   bombus_adult_cypermethrin_oral$LD50_data, bombus_adult_cypermethrin_contact$LD50_data,
                   bombus_larvae_acetamiprid_contact$LD50_data,
                   bombus_larvae_cypermethrin_contact$LD50_data)
ggtexttable(LD50_data, rows = NULL, theme = ttheme("lBlueWhite", base_size = 10))
