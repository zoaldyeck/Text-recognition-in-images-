
library(EBImage)
library(ForeCA)
#remove all the variables
rm(list=ls())

#k mean, unit=16, s=2,pNum=120000
k<-400
unit<-16
s<-4
pNum<-200000
all<-5000  #posters num, originals num


#grayresize read in all images in a file, grayscale and resize to be 96*96, then store
#the matrix into a list
GrayResize <- function (imgList) {
  matList <- list()
  imgNum <- min(length(imgList), all)
  
  for(i in 1:imgNum) {
    if(i%%1000==0){
      cat("Reading images",i,"\n")
    }
    #read in images, grayscale, resize to be 98*98, convert to matrix.store in matList
    magic <- readBin(imgList[i], what = 0L, n = 8, size = 1L, signed = FALSE)
    judgejpg <- isTRUE(all.equal(magic[1:2], c(0xFF, 0xD8)))
    judgepng <- isTRUE(all.equal(magic[1:8], c(0x89,0x50,0x4E,0x47,0x0D, 0x0A, 0x1A, 0x0A)))
    
    if(judgejpg || judgepng){
      img <- readImage(imgList[i])
      grayImg <- channel(img,"gray")
      grayImg <- resize(grayImg,96,96)
      gray_Mat <- as.matrix(grayImg)
      gray_Mat <- round(gray_Mat,digits=2)
     
      matList[[i]] <- gray_Mat
     
      
    }
  }
  return(matList)
}


#training+testing list
poster_list<-list()
origin_list<-list()
#poster.train.list, poster.test.list
poster_train_list<-list()
poster_test_list<-list()
#original.train.list,original.test.list
origin_train_list<-list()
origin_test_list<-list()
#poster+original train list, poster+original test list
train_list<-list()
test_list<-list()


#setwd("Path/to/posters")
setwd("/protestify_data/images/addlexpts/poster")
#poster_jpg, poster_png are string list of files names.
poster_jpg <- list.files(pattern="*.jpg")
poster_list <- GrayResize(poster_jpg)
#setwd("Path/to/originals")
setwd("/protestify_data/images/originals")
#origin_jpg, origin_png are string list of files names.
origin_jpg <- list.files(pattern = "*.jpg")
origin_list <- GrayResize(origin_jpg)



poster_num <- length(poster_list)
origin_num <- length(origin_list)



################################################################
###                 matrix of images list                    ###
################################################################
#randomly select 20% as test data
#divide the posters to test part(20%), and training+validation part(80%)
#randomly select 20% as test data
poster_test_index <- sample(1:poster_num, round(0.10*poster_num))
origin_test_index <- sample(1:origin_num, round(0.10*origin_num))

poster_test_list <- poster_list[poster_test_index]       #108
poster_train_list <- poster_list[-poster_test_index]     #432

origin_test_list <- origin_list[origin_test_index]       #74
origin_train_list <- origin_list[-origin_test_index]     #295

train_list <- c(poster_train_list,origin_train_list)  #727
test_list <- c(poster_test_list,origin_test_list)     #182
###############################################################
###    read images finished here       (cost 4min)          ###
###############################################################

###################################################################################
###    randomly select patches from all images, Use k-means to find k centers    ##
###################################################################################
#loop through all the unlabeled training images, randomly select 16*16 patches from 
#images, pNum(140000) patches in total are selected.
patch_extraction <- function(imgList,pNum,unit) {
  count <- 0
  patch_vec_list <- list()
  img_num <- length(imgList)
  
  # number of patch from each image
  num <- floor(pNum/img_num)
  
  for(i in 1:img_num) {
    for(j in 1:num) {
      #randomly generate row, col
      posRow <- sample(1:(96-unit),1)
      posCol <- sample(1:(96-unit),1)
      
      #unit=16, 16*16 patch->16*16 matrix->256 vector
      patch <- imgList[[i]][posRow:(posRow+unit-1),posCol:(posCol+unit-1)]
      
      patch_vec <- as.vector(patch)
      count <- count + 1
      
      patch_vec_list[[count]] <- patch_vec
    }
  }
  patch_mat <- matrix(unlist(patch_vec_list),nrow = unit*unit,byrow = FALSE)
  return(patch_mat)
}

#get the patches from train list
patch_Mat <- patch_extraction(train_list,pNum,unit)
###################################################################################
###  pNum=120000, select 120000 patches from 727 images,ends here take 3 min    ###
###################################################################################

##################################################################################
## normalize each patch vector by subtracting the means and then                ##
## dividing by the standard diviation                                           ##
##################################################################################
pre_process <- function(Mat) {
  #sd is the column variance
  sd <- sqrt(apply(Mat, 2, var))
  #var=0, replace it with 1
  #every vector - column mean
  sub_mean <- t(t(Mat) - colMeans(Mat))
  sd[sd==0]<-1
  normMat <- sweep(sub_mean,2,sd,FUN = "/")
  normMat<-round(normMat, digits=2)
  return(normMat)
}
normMat <- pre_process(patch_Mat)

############################################################
##                 ZCA whitening                          ##
############################################################
#reduce the correlation between vectors
ZCAwhite <- function(x,epsilon) {
  row <- nrow(x)
  col <- ncol(x)
  sigma <- x %*% t(x) / col
  duv <- svd(sigma)
  #  xRot <- t(duv$u) %*% x
  #  xTrilde <- t(duv$u[,1:col]) %*% x
  xPCA <- diag(1./sqrt((duv$d)+epsilon)) %*% t(duv$u) %*% x
  xZCA <- duv$u %*% xPCA
  xZCA<-round(xZCA,digits=2)
  return(xZCA)
}
matZCA <- ZCAwhite(normMat,0.01)


#norm_vec<-as.vector(normMat)
#PCA_vec<-princomp(norm_vec)
#matZCA<-matrix(PCA_vec$scores,nrow=nrow(normMat),byrow=FALSE)
############################################################
##          k-mean cluster k=800, take 15min              ##
############################################################
#pNum(120000) patches(vectors) in total, use k-means to find k(800) centers(vectors)
k_result <- kmeans(t(matZCA),k)
kcentroids <- k_result$centers
kcentroids <- t(kcentroids)



###########################################################################
##Uniformly select patches from images, convert to 256-dimensional       ##
## vectors, then convert to k-dimensional vectors, do locally sum        ##
## up to reduce the dimension, and get a 9k dimensional vector in the end##
###########################################################################
#matsplitter take in a 96*96 matrix(one image), every patch is 16*16, the leftmost and
#upmost of every adjacent patch has a distance of 4, one image->441 patches, every side
#was divided to 21 part((96-16)/4+1).
#since we need to do locally sum up later(combine 49 patches, 441 patches->9 patches), 
#we store the vectors in a different sequence
matsplitter<-function(M, u, s) {
  cv<-matrix(NA,nrow=unit*unit,ncol=441,byrow=FALSE)
  rn<-(96-u)/s+1     #21
  group<-floor(rn/3)  #7
  count<-0
  #one image is divided to 36 parts.
  for(r in 0:2){
    for(c in 0:2){
      #the row, col offset of every part
      r_off<-r*group*s+1
      c_off<-c*group*s+1
      #every part contains 49 pathces(overlap together), group=7
      for(x in 0:(group-1)){
        for(y in 0:(group-1)){
          #the row, col offset of every patch
          x_off<-r_off+x*s
          y_off<-c_off+y*s
          mat<-M[x_off:(x_off+u-1),y_off:(y_off+u-1)]
          vec<-as.vector(mat)
          count<-count+1
          #convert to vector, stores in cv
          cv[,count]<-vec
        }
      }
    }
  }
  return(cv)
} 
#calculate euclidean distance between 2 matrix: pMat, centroids 
#every vector in pMat is a patch, every vector in centroids is our k centers
#pMat(256*1764), centroids(256*800)
euclidean_dist<-function(pMat,centroids){
  G<-t(pMat)%*%centroids
  Ss<-diag(t(pMat)%*%pMat)
  Rr<-diag(t(centroids)%*%centroids)
  n<-ncol(pMat)
  S<-matrix(Ss,nrow=n,ncol=k,byrow=FALSE)
  R<-matrix(Rr,nrow=n,ncol=k,byrow=TRUE)
  
  dis<-sqrt(S+R-2*G)
  return(dis)
}
#img_feature take in one image, and k centroids, split into 1764 patches, normalize 
#calculate the euclidean distance with the centroids
feature_Generator <- function(inputList,centroids) {
  num <- length(inputList)
  fMat <- matrix(NA,nrow=9*k,ncol=num,byrow=FALSE)
  count<-0
  for(i in 1:num) {
    if(i%%100==0){
      cat("Converting images",i,"\n")
    }
    
    temp<-inputList[[i]]
    number<-ncol(imageData(temp))
    
    if(!is.null(number)){
    pMat <- matsplitter(temp,unit,s)
    
    pMat <- pre_process(pMat)
    # pMat <- ZCAwhite(pMat,0.01)
    dis<-euclidean_dist(pMat,centroids)
    #udis is the average of the distance between one vector(in image) with k centroids
    udis<-rowMeans(dis)
    dif_dis<-dis-udis
    #max(0,dis-mean)
    dif_dis[dif_dis<0]<-0
    #it's called k-means(triangle) clustering
    #combine 196 patches together and get an average
    fvList <- matrix(,nrow=k,ncol=9,byrow=FALSE)
    for(j in 0:8){
      row<-j*49+1
      temp<-dif_dis[row:(row+48),]
      local_sum<-colMeans(temp)
      fvList[,(j+1)] <- local_sum
    }
    count<-count+1
    vec<-as.vector(fvList)
    fMat[,count]<-vec
    }
  }
  fMat<-fMat[,1:count]
  return(fMat)
}

poster_train__feature_mat <- feature_Generator(poster_train_list,kcentroids)
origin_train_feature_mat <- feature_Generator(origin_train_list,kcentroids)

poster_test__feature_mat <- feature_Generator(poster_test_list,kcentroids)
origin_test_feature_mat <- feature_Generator(origin_test_list,kcentroids)

###################################################################################
#### all data transition ends here, takes 10 min                              #####
###################################################################################
#add label:poster=1. original=-1
poster_train <- rbind(poster_train__feature_mat,+1)
origin_train <- rbind(origin_train_feature_mat,-1)
poster_test <- rbind(poster_test__feature_mat,+1)
origin_test <- rbind(origin_test_feature_mat,-1)

train <- cbind(poster_train,origin_train)
test <- cbind(poster_test,origin_test)





####################################################################################
####################     5-fold on training data, average of 81.6%             #######
####################################################################################
#This part is for 5-fold cross validation, all the fk function learning is based on
#training data, divide training data to 5 part, train 5 times on 80%, validate on the rest.
#shuffle before validation.
train_num<-ncol(train)
shuffle<-sample(1:train_num, train_num)
shuffle_train<-train[,shuffle]

fold<-5
fen<-floor(train_num/fold)
for(i in 1:fold){
  #get 20% out to be validation
  val_index<-c(((i-1)*fen+1):(i*fen))
  
  val_sample<-shuffle_train[,val_index]
  #divide to val_x and val_y
  val_x<-val_sample[1:(9*k),]
  val_x<-t(val_x)
  val_y<-val_sample[(9*k+1),]
  #get the rest 80% to be training data
  train_sample<-shuffle_train[,-val_index]
  #divide to train_x and train_y
  train_x<-train_sample[1:(9*k),]
  train_x<-t(train_x)
  train_y<-train_sample[(9*k+1),]
  
  library(kernlab)
  #use svm on training data to get a model
  model <- ksvm(train_x,train_y,type="C-svc",kernel='vanilladot',C=100,scaled=c())
  #use the model to predict y_test_pred from val_x, then compare to val_y
   y_test_pred<-predict(model,val_x)
  accuracy<-sum(y_test_pred==val_y)/length(y_test_pred)
  cat("K=1600,unit=16,s=4,accuracy=",accuracy)
}





###################################################################################
############           test        88.2%                                      ###
###################################################################################
#do the same thing on test data
test_x<-test[1:(9*k),]
test_x<-t(test_x)
test_y<-test[(9*k+1),]

y_test_pred<-predict(model,test_x)

test_accuracy<-sum(y_test_pred==test_y)/length(y_test_pred)
cat("testing:K=1600,unit=16,s=2,accuracy=",test_accuracy)





###########################################################################
##           test groups 85%~91% for originals, 70%~83% for posters      ##
###########################################################################
#This part is for testing images in one file by setting the path, and read in all
#the images.
#setwd("E:/")
#setwd("E:/poster/")
#group_jpg, group_png are string list of files names.
#group_jpg <- list.files(pattern="*.jpg")
#group_png <- list.files(pattern="*.png")
#group_jpg_list, group_png_list store matrixs of jpg and png images.
#group_jpg_list <- list()
#group_png_list <- list()
#group_jpg_list <- GrayResize(group_jpg)
#group_png_list <- GrayResize(group_png)
#group_list <- c(group_jpg_list,group_png_list)
#group_test_mat <- feature_Generator(group_list,kcentroids)
#group_test_pred<-predict(model,t(group_test_mat))
#sum(group_test_pred+1)/(2*length(group_test_pred))   #count the accuracy for the originals
#1+sum(group_test_pred-1)/(2*length(group_test_pred))   #count the accuracy for the posters


#setwd("E:/origin/")
#setwd("E:/poster/poster/")
#group_jpg, group_png are string list of files names.
#group_jpg <- list.files(pattern="*.jpg")
#group_png <- list.files(pattern="*.png")
#group_jpg_list, group_png_list store matrixs of jpg and png images.
#group_jpg_list <- list()
#group_png_list <- list()
#group_jpg_list <- GrayResize(group_jpg)
#group_png_list <- GrayResize(group_png)
#group_list <- c(group_jpg_list,group_png_list)
#group_test_mat <- feature_Generator(group_list,kcentroids)
#group_test_pred<-predict(model,t(group_test_mat))
#1-sum(group_test_pred+1)/(2*length(group_test_pred))   #count the accuracy for the originals
#sum(group_test_pred-1)/(2*length(group_test_pred))   #count the accuracy for the posters
