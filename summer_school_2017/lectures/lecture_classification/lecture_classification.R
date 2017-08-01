# Lecture "Classification". R script.

library(ggplot2)
library(class)
library(foreign)
library(pROC)

data("iris")
summary(iris)

# Glance look

ggplot(iris, aes(x=Sepal.Length, Petal.Length, color=Species)) +
  geom_point() + theme_bw()

# Splitting to train/valid

set.seed(1)
training <- sample(1:150, 120)
validation <- (1:150)[-training]

train <- iris[training, ]
valid <- iris[validation, ]

ggplot(train, aes(x=Sepal.Length, Petal.Length, color=Species)) +
  geom_point() + theme_bw() +
  geom_point(data=valid, mapping = aes(color=NA))

# K-nn

prediction <- knn(train[, 1:4], valid[, 1:4], train$Species, k=5)
prediction
table(valid$Species, prediction, dnn=c("actual", "predicted"))


# Diabetes data

data <- read.arff("diabetes.arff.txt")
head(data)

# Splitting to train/valid 

set.seed(1)
validation <- sample(1:nrow(data), nrow(data) %/% 5)
training <- (1:nrow(data))[-validation]

train <- data[training, ]
valid <- data[validation, ]

# Generalised linear regression for logistic regression purpose

model <- glm(class ~ ., 
             family = binomial, data=train)
plot(model)
prediction <- predict(model, valid, type="response")
curve(prediction)

prediction <- ifelse(prediction > 0.243578190370293, "tested_positive", "tested_negative")
table(valid$class, prediction)

# ROC curve

probs <- predict(model, valid, type=c("response"))
valid$probs <- probs
g <- roc(class ~ probs, data = valid)
plot(g)


for (i in 1:154) {
  print(paste(
    g$sensitivities[i],
    g$specificities[i],
    g$thresholds[i],
    collapse = " "
  ))
}
