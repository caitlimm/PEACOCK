library(mlr)
# sampledata is from MCF7 and contains standardized variables ending in Z, which were standardized based on the mean 
##### and SD values from all MCF7 features values for enhancer-gene links in cislinks file (~17M enhancer-gene pairs)

#final model
vars<-variable.names(sampledata)[-1]
vars<-vars[-1]
test<-sampledata[vars]
test$truelink=7 #this is a nonsense value just to construct the matrix
ytest=matrix(test$truelink)
xtest=model.matrix(test$truelink ~., data=test)[,-1]
ridge.pred=predict(final, newx=xtest, type="response")
sampledata$score<-ridge.pred #each value in the sample data now has a score based on the final model

#alternate final model
rf_predict = predict(final_alternate, type="prob", newdata=sampledata, probability=TRUE)
sampledata$score<-rf_predict$data$prob.1 #each value in the sample data now has a score based on the final model