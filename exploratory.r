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
library(Metrics)
library(randomForestSRC)

#Distribution of follow up time
brca <- fread("brca.csv", na.strings = "")
brcaclinic <- brca[, 2:9]
brca_prot <- brca[, 16:213]
brca_df <- cbind(brcaclinic, brca_prot %>% select_if(~ !any(is.na(.))))
brca_df$race[is.na(brca_df$race)] <- "Unknown"
brca_df$subtype[is.na(brca_df$subtype)] <- "Unknown"
brca_df <- brca_df[brca_df$race!="AMERICAN INDIAN OR ALASKA NATIVE",] 
brca_df$ajcc_pathologic_tumor_stage[is.na(brca_df$ajcc_pathologic_tumor_stage)] <- "Unknown"
brca_data <- brca_df %>% dummy_cols()
brca_data <- brca_data %>% select_if(~ !is.character(.)) %>% clean_names()

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

brca_boost_train <- train  

set.seed(123)
cox_boost <- 
  boosting_core( 
  Surv(os_time, os)~.,
  rate = 0.1, 
  num_iter = 3000,
  control_method = "BIC",
  data = brca_boost_train)

boost_df <- 
  train %>% 
  select(c(names(cox_boost$coefficients[cox_boost$coefficients > 0]), os_time, os)) %>%
  select( -prex1, -eif4e, -ecadherin, -cmyc, -claudin7, -ar)

cox_fit_train <- coxph(Surv(os_time, os)~. , data = boost_df, x=T)
testPH <- cox.zph(cox_fit_train)
testPH

cox_fit_train %>% summary()

train_lp <- cox_fit_train$linear.predictors
train_lp_pred <- predict(cox_fit_train, test, type = "lp")
rmse(train_lp, train_lp_pred)


cox_fit_test <- coxph(Surv(os_time, os)~., data = test %>% select(names(boost_df)), x = T)
cox_fit_test %>% summary()

#31.1
set.seed(123)
rf_fit_train <- rfsrc(Surv(os_time, os)~. , data = boost_df, ntree = 1000, mtry = 4)
rf_fit_train 


#30.59
set.seed(123)
rf_fit_train <- rfsrc(Surv(os_time, os)~. , data = boost_df, ntree = 2000, mtry = 5)
rf_fit_train 
cindex(rf_fit_train)

#30.85
set.seed(123)
rf_fit_train <- rfsrc(Surv(os_time, os)~. , data = boost_df, ntree = 3000, mtry = 6)
rf_fit_train 

set.seed(123)
rf_fit_test <- rfsrc(Surv(os_time, os)~. , data = test %>% select(names(boost_df)), ntree = 2000, mtry = 5)
rf_fit_test 
cindex(rf_fit_test)

