# Lecture "Regression". R script

library(ggplot2)

# Quick data overview
ageData <- read.csv("age_data.tsv", sep="\t")
head(ageData)

ggplot(data=ageData, aes(x=Age.End, y=Male.Cases)) +
  geom_point(size=3) + theme_bw()


# Splitting data into train/validation
set.seed(1)
training <- sample(1:nrow(ageData), 12)
validation <- (1:nrow(ageData))[-training]

ageData$set <- NA
ageData[training, "set"] <- "training"
ageData[validation, "set"] <- "validation"

head(ageData[ , c(1:3, 7)])

ggplot(data=ageData, aes(x=Age.End, y=Male.Cases, color=set)) +
  geom_point(size=3) + theme_bw()

# linear regression

ggplot(data=ageData, aes(x=Age.End, y=Male.Cases, color=set)) +
  geom_point(size=3) + theme_bw() + 
  geom_smooth(data=subset(ageData, set == "training"),
              method="lm", formula = y ~ x, se=F, fullrange=T)

# visualising error

train <- subset(ageData, set == "training")

with(train, {
  fit <- lm(Male.Cases ~ Age.End)
  prediction <- predict(fit, data.frame('Age.End' = ageData$Age.End))
  predicted <- cbind(ageData, prediction)
  
  ggplot(data=predicted, aes(x=Age.End, y=Male.Cases, color=set)) +
    geom_point(size=3) + 
    geom_point(aes(y = prediction), shape = 1, size=3) +
    geom_segment(aes(xend = Age.End, yend = prediction), alpha = .8) +
    geom_smooth(data=subset(ageData, set == "training"),
                method="lm", formula = y ~ x, se=F, fullrange=T) +
    theme_bw()
  
})

## Polynomial regression

degree <- 2
ggplot(data=ageData, aes(x=Age.End, y=Male.Cases, color=set)) +
  geom_point(size=3) + theme_bw() + 
  geom_smooth(data=subset(ageData, set == "training"),
              method="lm", formula = y ~ poly(x, degree), se=F, fullrange=T)

## Overfitting

with(train, {
  results <- lapply(1:10, function(degree) {
    fit <- lm(Male.Cases ~ poly(Age.End, degree=degree, raw=T))
    prediction <- predict(fit, data.frame('Age.End' = ageData$Age.End))
    se <- (ageData$Male.Cases - prediction)^2
    
    rmseTrain <- sqrt(mean(se[training]))
    rmseValid <- sqrt(mean(se[validation]))
    return(c(degree, rmseTrain, rmseValid))
  })
  
  results <- do.call(rbind, results)
  colnames(results) <- c("Degree", "SSE Train", "SSE Validation")
  toPlot <- rbind(
    data.frame(degree=results[, "Degree"], SSE=results[, "SSE Train"], dataset="Train"),
    data.frame(degree=results[, "Degree"], SSE=results[, "SSE Validation"], dataset="Validation")
  )
  
  ggplot(data=toPlot, aes(x=degree, y=SSE, color=dataset)) +
    geom_point(size=3) + geom_line(size=2) + ggtitle("RMSE Plot") + theme_bw()
  
})