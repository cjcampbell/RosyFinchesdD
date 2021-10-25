library(MuMIn)
library(car)
df <- myResults %>% 
  filter(Season == "Winter") %>% 
  filter(!is.na(Sex), !is.na(Size), !is.na(Region_long))

m_global <- lm(probDistTraveledAtThreshold_0.25~Sex+Size+Region_long, data = df)
summary(m_global) 
car::vif(m_global)

options(na.action = "na.fail")
m_dredged <- dredge(m_global)
# View(m_dredged)

# Model weight (after discarding noncompeting models)
m_dredged_clean <- subset(m_dredged, delta==0|delta>2)
m_dredged_clean

# Check out top model
topModel <- lm(probDistTraveledAtThreshold_0.25~Size+Region_long, data = df)
# How much variance is explained?
out <- performance::r2(topModel) 


# Put model estimates into temporary data.frames:
model1Frame <- data.frame(Variable = rownames(summary(topModel)$coef),
                          Coefficient = summary(topModel)$coef[, 1],
                          SE = summary(topModel)$coef[, 2],
                          modelName = "Top Model", row.names = NULL)
model2Frame <- data.frame(Variable = rownames(summary(m_global)$coef),
                          Coefficient = summary(m_global)$coef[, 1],
                          SE = summary(m_global)$coef[, 2],
                          modelName = "Global Model",  row.names = NULL)

# Combine these data.frames
allModeldf <- data.frame(rbind(model1Frame, model2Frame)) %>% 
  dplyr::mutate(
    modelName = factor(modelName, levels = c("Global Model","Top Model")),
    Variable = case_when(
      Variable == "Region_longNorth-central" ~ "Hib. Region (N-C)",
      Variable == "SizeSmall" ~ "Hib. Size (Small)",
      Variable == "SexMale" ~ "Sex (Male)",
      TRUE ~ Variable
      ),
      Variable = factor(Variable, levels = rev(c("Hib. Region (N-C)", "Hib. Size (Small)", "Sex (Male)", "(Intercept)")))
    )

save(m_global, topModel, allModeldf,  file = file.path(wd$bin, "modelResults.Rdata"))
