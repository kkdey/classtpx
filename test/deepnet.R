
####  Deep Learning with the mxnet package ################################

data(Sonar, package="mlbench")
Sonar[,61] = as.numeric(Sonar[,61])-1
train.ind = c(1:50, 100:150)
train.x = data.matrix(Sonar[train.ind, 1:60])
train.y = Sonar[train.ind, 61]
test.x = data.matrix(Sonar[-train.ind, 1:60])
test.y = Sonar[-train.ind, 61]


mx.set.seed(0)
model <- mx.mlp(train.x, train.y, hidden_node=50, out_node=2, out_activation="softmax",
                num.round=50, array.batch.size=15, learning.rate=0.1, momentum=0.8, 
                eval.metric=mx.metric.accuracy)

graph.viz(model$symbol$as.json())
preds = predict(model, test.x)

pred.label = max.col(t(preds))-1
table(pred.label, test.y)

training.data.frame <- data.frame(cbind(train.y, train.x));
model <- svm(as.factor(train.y) ~ ., data = training.data.frame)

test_class_svm <- predict(model, test.x)
tab_class_svm <- table(test_class_svm, test.y)

data <- mx.symbol.Variable("data")
fc1 <- mx.symbol.FullyConnected(data, name="fc1", num_hidden=128)
act1 <- mx.symbol.Activation(fc1, name="relu1", act_type="relu")
fc2 <- mx.symbol.FullyConnected(act1, name="fc2", num_hidden=64)
act2 <- mx.symbol.Activation(fc2, name="relu2", act_type="relu")
fc3 <- mx.symbol.FullyConnected(act2, name="fc3", num_hidden=10)
softmax <- mx.symbol.SoftmaxOutput(fc3, name="sm")


data <- mx.symbol.Variable("data")
fc1 <- mx.symbol.FullyConnected(data, name="fc1", num_hidden=128)
act1 <- mx.symbol.Activation(fc1, name="tanh1", act_type="tanh")
fc2 <- mx.symbol.FullyConnected(act1, name="fc2", num_hidden=64)
act2 <- mx.symbol.Activation(fc2, name="tanh2", act_type="tanh")
fc3 <- mx.symbol.FullyConnected(act2, name="fc3", num_hidden=32)
act3 <- mx.symbol.Activation(fc3, name="tanh3", act_type="tanh")
fc4 <- mx.symbol.FullyConnected(act3, name="fc4", num_hidden=16)
act4 <- mx.symbol.Activation(fc4, name="tanh4", act_type="tanh")
fc5 <- mx.symbol.FullyConnected(act4, name="fc5", num_hidden=8)
act5 <- mx.symbol.Activation(fc5, name="tanh5", act_type="tanh")
fc6 <- mx.symbol.FullyConnected(act5, name="fc6", num_hidden=2)
softmax <- mx.symbol.SoftmaxOutput(fc6, name="sm")

devices <- mx.cpu()

mx.set.seed(0)
model <- mx.model.FeedForward.create(softmax, X=train.x, y=train.y,
                                     ctx=devices, num.round=10, array.batch.size=10,
                                     learning.rate=0.07, momentum=0.9,  eval.metric=mx.metric.accuracy,
                                     initializer=mx.init.uniform(0.07),
                                     epoch.end.callback=mx.callback.log.train.metric(100))

graph.viz(model$symbol$as.json())




library(mxnet)
data(BostonHousing, package="mlbench")
train.ind = seq(1, 506, 3)
train.x = data.matrix(BostonHousing[train.ind, -14])
train.y = BostonHousing[train.ind, 14]
test.x = data.matrix(BostonHousing[-train.ind, -14])
test.y = BostonHousing[-train.ind, 14]

data <- mx.symbol.Variable("data")

fc1 <- mx.symbol.FullyConnected(data, num_hidden=2)

lro <- mx.symbol.LinearRegressionOutput(fc1)
mx.set.seed(0)
model <- mx.model.FeedForward.create(lro, X=train.x, y=as.array(train.y),
                                     ctx=mx.cpu(), num.round=50, array.batch.size=20,
                                     learning.rate=2e-6, momentum=0.9, eval.metric=mx.metric.rmse)



########  Deep Learning with the h2o package ############################

library(h2o)
localH2O <- h2o.init(ip = "localhost", port = 54321, startH2O = TRUE)
data(Sonar, package="mlbench")
train.ind = c(1:50, 100:150)
train.x = data.matrix(Sonar[train.ind, 1:60])
train.y = Sonar[train.ind, 61]

training_dat <- data.frame(cbind(train.y, train.x));


h2o.init()
prosPath <- system.file("extdata", "prostate.csv", package="h2o")
prostate.hex <- h2o.uploadFile(path = prosPath)
as.data.frame(prostate.hex)

h2o.dat <- as.h2o(training_dat, destination_frame = "dat")

model <- 
  h2o.deeplearning(x = 2:60,  # column numbers for predictors
                   y = 1,   # column number for label
                   training_frame = as.h2o(training_dat, destination_frame = "dat"),
                  # data = dat, # data in H2O format
                   activation = "TanhWithDropout", # or 'Tanh'
                   input_dropout_ratio = 0.2, # % of inputs dropout
                   hidden_dropout_ratios = c(0.5,0.5,0.5), # % for nodes dropout
                   balance_classes = TRUE, 
                   hidden = c(50,50,50), # three layers of 50 nodes
                   epochs = 100) # max. no. of epochs