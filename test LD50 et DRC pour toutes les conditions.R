# Herluin Saillard, with Antoine Gekière (Laboratory of Zoology, University of Mons, Belgium) - April 2024
library(readxl)
library(dplyr)
library(drc)
library(ggplot2)
library(ggpubr)
library(patchwork)

# Load the data ----

data <- bombus_acute_et_oral_test


#-------------------------   Dose Response Curve and LC50 dataframe   ----
# OECD (2006), Current approaches in the statistical analysis of ecotoxicity data: a guidance to application
# Ritz et al (2015), Dose-response analysis using R

analyze <- function(data, species_var, stage_var, pesticide_var, exposition_var) {
  
## Antoine Gekière's code DRC - April 2024 ----
  
  # Filter the data for the specific substance and remove NA
  mortality_data <- data %>%
    filter((species = species_var & stage = stage_var & pesticide == pesticide_var | pesticide == "Negative") & exposition == exposition_var) %>% # need different names (eg. "stage = stage" don't work)
    filter(!is.na(censored))
  
  # Generate dataset for analysis
  data_drc <- mortality_data %>%
    group_by(conc_nomil) %>%
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
  #              (µg / 40 µL for oral)
  #              (µg / 02 µL for contact)
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
    labs(title = paste("DRC for bombus /", stage_var, "/", exposition_var, "/", pesticide_var),
         x = "Concentration (mg / L)",
         y = "Mortality") +
    theme_bw()
  
## Create a dataframe with LC50, t-value and p-value ----
  LC50_data <- subset(data.frame
                      (Parameter = rownames(sum_drc$coefficients),
                        Estimate_LC50 = sum_drc$coefficients[, "Estimate"],
                        Standard_Error = sum_drc$coefficients[, "Std. Error"],
                        t_value = sum_drc$coefficients[, "t-value"],
                        p_value = sum_drc$coefficients[, "p-value"]),
                      Parameter != "Slope:(Intercept)")
  LC50_data$Parameter <- paste("Bombus", stage_var, pesticide_var, exposition_var)  # replace name of the parameter by "LC 50 [substance name]"
  
## Var lists ----
  return(list(drc = drc, plot = plot, LD50 = LD50, standard_Error = standard_error, sum_drc = sum_drc, LC50_data = LC50_data))
}

#-------------------------   Modalities ----

bombus_adult_acetamiprid_oral <- analyze(data, "bombus", "adult", "Acetamiprid", "oral")
bombus_larvae_acetamiprid_oral <- analyze(data, "bombus", "larvae", "Acetamiprid", "oral")

bombus_adult_cypermethrin_oral <- analyze(data, "bombus", "adult", "Cypermethrin", "oral")
bombus_larvae_cypermethrin_oral <- analyze(data, "bombus", "larvae", "Cypermethrin", "oral")

bombus_adult_cypermethrin_contact <- analyze(data, "bombus", "adult", "Cypermethrin", "contact")
bombus_larvae_cypermethrin_contact <- analyze(data, "bombus", "larvae", "Cypermethrin", "contact")


#-------------------------   Display drc plot and text table     ----

# Display the dose response curve
print(bombus_adult_acetamiprid_oral$plot)
print(bombus_larvae_acetamiprid_oral$plot)

print(bombus_adult_cypermethrin_oral$plot)
print(bombus_larvae_cypermethrin_oral$plot)

print(bombus_adult_cypermethrin_contact$plot)
print(bombus_larvae_cypermethrin_contact$plot)

# Create and display the text table
LC50_data <- rbind(bombus_adult_acetamiprid_oral$LC50_data, bombus_larvae_acetamiprid_oral$LC50_data,
                   bombus_adult_cypermethrin_oral$LC50_data, bombus_larvae_cypermethrin_oral$LC50_data,
                   bombus_adult_cypermethrin_contact$LC50_data, bombus_larvae_cypermethrin_contact$LC50_data)
ggtexttable(LC50_data, rows = NULL, theme = ttheme("lBlueWhite", base_size = 10))


#-------------------------   LD50 bar plot with all treatments     ----

# Create a dataframe with LD 50 and standard error
# Use Copilot for auto write all modalities
LD50_bombus_adult_acetamiprid_oral <- bombus_adult_acetamiprid_oral$LD50
standard_error_bombus_adult_acetamiprid_oral <- bombus_adult_acetamiprid_oral$standard_Error
LD50_bombus_larvae_acetamiprid_oral <- bombus_larvae_acetamiprid_oral$LD50
standard_error_bombus_larvae_acetamiprid_oral <- bombus_larvae_acetamiprid_oral$standard_Error

LD50_bombus_adult_cypermethrin_oral <- bombus_adult_cypermethrin_oral$LD50
standard_error_bombus_adult_cypermethrin_oral <- bombus_adult_cypermethrin_oral$standard_Error
LD50_bombus_larvae_cypermethrin_oral <- bombus_larvae_cypermethrin_oral$LD50
standard_error_bombus_larvae_cypermethrin_oral <- bombus_larvae_cypermethrin_oral$standard_Error

LD50_bombus_adult_cypermethrin_contact <- bombus_adult_cypermethrin_contact$LD50
standard_error_bombus_adult_cypermethrin_contact <- bombus_adult_cypermethrin_contact$standard_Error
LD50_bombus_larvae_cypermethrin_contact <- bombus_larvae_cypermethrin_contact$LD50
standard_error_bombus_larvae_cypermethrin_contact <- bombus_larvae_cypermethrin_contact$standard_Error

LD50_data <-
  data.frame(
    modality = c("Bombus adult acetamiprid oral", "Bombus larvae acetamiprid oral",
                 "Bombus adult cypermethrin oral", "Bombus larvae cypermethrin oral",
                 "Bombus adult cypermethrin contact", "Bombus larvae cypermethrin contact"),
    LD50 = c(LD50_bombus_adult_acetamiprid_oral, LD50_bombus_larvae_acetamiprid_oral,
             LD50_bombus_adult_cypermethrin_oral, LD50_bombus_larvae_cypermethrin_oral,
             LD50_bombus_adult_cypermethrin_contact, LD50_bombus_larvae_cypermethrin_contact),
    standard_error = c(standard_error_bombus_adult_acetamiprid_oral, standard_error_bombus_larvae_acetamiprid_oral,
                       standard_error_bombus_adult_cypermethrin_oral, standard_error_bombus_larvae_cypermethrin_oral,
                       standard_error_bombus_adult_cypermethrin_contact, standard_error_bombus_larvae_cypermethrin_contact))

# Create the bar plot

bar_plot <-
  ggplot(LD50_data, aes(x = modality, y = LD50)) +
  geom_col(aes(fill = modality),
           color = "black") +
  scale_fill_manual(values = c("Bombus adult acetamiprid oral" = "#FA528D", "Bombus larvae acetamiprid oral" = "#F5B041",
                               "Bombus adult cypermethrin oral" = "#E74C3C", "Bombus larvae cypermethrin oral" = "#F4D03F",
                               "Bombus adult cypermethrin contact" = "#C0392B", "Bombus larvae cypermethrin contact" = "#F1C40F")) + #choice colors for all modalities
  geom_errorbar(aes(ymin = LD50 - standard_error, ymax = LD50 + standard_error),
                width = 0.2, position = position_dodge(0.9)) +
  labs(title = "LD 50 for Bombus Terrestris", x = "Modality", y = "LD 50 (µg / bee)") +  # Add axis titles
  theme_bw()+
  theme(panel.grid = element_line(color = "grey"))
print(bar_plot)
