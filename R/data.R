#' Simulated case control dataset for 6856 stroke cases and 6856 stroke controls
#'
#'
#' Dataset containing simulated data on risk factors for 6856 stroke cases and 6856 stroke control, based on risk factors and associations in the INTERSTROKE study
#'
#'
#' @format A data frame with 13712 rows and 19 variables:
#' \describe{
#'   \item{region}{Geographic region, 1: Western Europe, 2: Eastern/central Europe/Middle East 3: Africa, 4: South Asia, 5: China, 6: South East Asia, 7: South America}
#'   \item{case}{case control status, (1 for stroke cases)}
#'   \item{sex}{Gender of individual, 0: male, 1:female}
#'   \item{age}{Age of individual}
#'   \item{smoking}{Smoking status, 0: Never, 1: Current}
#'   \item{stress}{1: sometimes stressed, 0: never stressed}
#'   \item{waist_hip_ratio}{Waist hip ratio}
#'   \item{exercise}{Physical Activity.  1: mainly inactive, 0: mainly active}
#'   \item{alcohol}{Alcohol history and frequency, 1:never, 2:low/moderate, 3:high intake}
#'   \item{diabetes}{Diabetes, 0: No, 1: Yes}
#'   \item{diet}{Healthy eating score (higher is better)}
#'   \item{early_stage_heart_disease}{presence of risk factors for heart disease. 0: No, 1: yes}
#'   \item{lipids}{Ratio of Apolipoprotein B to Apolipoprotein A}
#'   \item{education}{Years of education.  1: No education, 2: 1-8 years, 3:9-12 years, 3:Technical college, 4: University}
#'   \item{high_blood_pressure}{Diagnosed hypertension: 0 No, 1: yes}
#'   \item{weights}{weights that are proportional to inverse sampling probabilities.  We have scaled the weights to be 0.0035 for a case and 0.9965 for a control to reflect any approximate incidence of 1 serious stroke in every 0.9965/0.0035 person years in the population}
#'   \item{time}{simulated time variable (for illustrating survival models)}
#'   \item{event}{simulated event indicator (0 if censored, 1 if event happened): for illustrating survival models}
#'   \item{strata}{Strata number based on sex and region.   For illustrating conditional regression}
#'
#'
#' }
#' @source Data simulated based on relationships described in \url{https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(16)30506-2/fulltext}
"stroke_reduced"



#' Simulated case control dataset for 5000 cases (individuals with chronic cough) and 5000 controls
#' @format A data frame with 10000 rows and 4 variables:
#' \describe{
#'   \item{y}{Chronic Cough, 1: Yes, 0: No}
#'   \item{urban.rural}{1: resident in urban setting, 0: resident in rural setting}
#'   \item{smoking.category}{Smoking level: 1 No smoker, 2: ex smoker, 3: 1-9 cigarettes per day, 4: 10-19 cigarettes per day, 4:>= 20 cigarettes per day}
#'   \item{occupational.exposure}{Exposed to dust/gas at work.  1: Yes, 0: no}
#'
#'   }
#' @source Data simulated based on "Sequential and average attributable fractions as aids in the selection of preventive strategies." Journal of clinical epidemiology 48, no. 5 (1995): 645-655.
"Hordaland_data"
