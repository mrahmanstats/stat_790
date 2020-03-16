library(data.table)
library(janitor)
library(tidyverse)
library(Amelia)
library(survival)
library(gtsummary)
library(survminer)
library(corrplot)
library(pec)
library(SurvBoost)
library(riskRegression)
library(fastDummies)

#Distribution of follow up time
brca <- fread("brca.csv", na.strings = "")
brcaclinic <- brca[, 2:9]
brca_prot <- brca[, 16:213]
brca_df <- cbind(brcaclinic, brca_prot %>% select_if(~ !any(is.na(.))))
brca_df$race[is.na(brca_df$race)] <- "Unknown"
brca_df$subtype[is.na(brca_df$subtype)] <- "Unknown"
brca_df <- brca_df[brca_df$race!="AMERICAN INDIAN OR ALASKA NATIVE",] 
brca_df$ajcc_pathologic_tumor_stage[is.na(brca_df$ajcc_pathologic_tumor_stage)] <- "Unknown"
brca_data <- brca_df 
# %>% dummy_cols()
# brca_data <- brca_data %>% select_if(~ !is.character(.))

ggplot(data = brca, aes(x = OS.time)) +
  geom_histogram(aes(fill = factor(OS))) +
  theme_classic() +
  labs(fill = "Status")

brca <- brca %>% rename("tumor_stage" = ajcc_pathologic_tumor_stage)

brca_summ <- brca[,1:15]
brca_summ <- brca_summ %>% dplyr::select(-patient, -PFI, -PFI.time, -DSS, -DSS.time, -DFI, -DFI.time)
tbl_summary(brca_summ, by = OS) %>% add_p() %>% bold_labels() %>% italicize_levels() %>% bold_p()

ggsurvplot(
  fit = survfit(Surv(OS.time, OS) ~ 1, data = brca), 
  xlab = "Days", 
  ylab = "Overall survival probability",
  pval = T,
  risk.table = T,
  surv.median.line = "hv")

ggsurvplot(
  fit = survfit(Surv(OS.time, OS) ~ tumor_stage, data = brca), 
  xlab = "Days", 
  ylab = "Overall survival probability",
  pval = T,
  risk.table = T)

# Random sample indexes
set.seed(123)
train_index <- sample(1:nrow(brca_data), 0.8 * nrow(brca_data))
test_index <- setdiff(1:nrow(brca_data), train_index)

train <- brca_data[train_index,]
test <- brca_data[test_index,]


cox_fit_1 <- coxph(Surv(OS.time, OS)~., data = train, x=T)
testPH <- cox.zph(cox_fit_1)
testPH

brca_boost <- train %>% dplyr::select(-ajcc_pathologic_tumor_stage, -histological_type, -race, -subtype, -age, -gender)

cox_boost <- boosting_core( 
  Surv(OS.time, OS)~.,
  rate = 0.2, 
  data = brca_boost)

boost_df <- train %>% dplyr::select(c(names(cox_boost$coefficients[cox_boost$coefficients > 0]), OS.time, OS))

cox_fit_2 <- coxph(Surv(OS.time, OS)~.-ECADHERIN, data = boost_df, x=T)
testPH <- cox.zph(cox_fit_2)
testPH

cox_fit_2 %>% summary()

cox_fit_2 <- coxph(Surv(OS.time, OS)~. -ajcc_pathologic_tumor_stage -subtype, data = train[,1:8], x=T)
testPH <- cox.zph(cox_fit_2)
testPH

cindex(cox_fit_1, Surv(OS.time, OS)~., data = brca_df[,1:8])
cindex(cox_fit_2, Surv(OS.time, OS)~ age + gender + race, data = brca_df[,1:8])

cox_fit2 <- coxph(Surv(OS.time, OS)~. -subtype,
                  data = brca_df)








