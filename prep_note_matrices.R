library(RColorBrewer)
library(tidyverse)
library(ggpmisc)
library(dplyr)
library(cowplot)
library(vegan)
library(ggpubr)
library(emmeans)
library(nlme)
library(car)
library(clusrank)

# whole song measurements for all songs in the dataset
d <- read_csv("all_song_data.csv",col_names=TRUE)

ids <- c("41-13-1","49-12-1","53-2-1","54-2-2","55-4-3","55-5-2","55-6-1","56-1-3","56-3-2","63-1-1")
treatment <- c("sham","ablation","ablation","sham","ablation","ablation","sham","sham","ablation","sham")

# read in all note matrices from the pre surgery recording period - each animal will have its own list
for(i in 1:length(ids)){
  list_pre <- list.files(path="~/steg_muscle_ablation/note_matrices/pre-surgery/",
                         pattern=paste0(ids[i],"*"),
                         full.names=TRUE)
  assign(paste0("pre_nm_",ids[i]),
         lapply(list_pre,read.csv))
}

# read in all note matrices from the post surgery recording period - each animal will have its own list
for(i in 1:length(ids)){
  list_post <- list.files(path="~/steg_muscle_ablation/note_matrices/post-surgery/",
                          pattern=paste0(ids[i],"*"),
                          full.names=TRUE)
  assign(paste0("post_nm_",ids[i]),
         lapply(list_post,read.csv))
}

# add frequency bandwidth as a variable to all lists
n <- c("`pre_nm_41-13-1`","`pre_nm_49-12-1`","`pre_nm_53-2-1`","`pre_nm_54-2-2`","`pre_nm_55-4-3`","`pre_nm_55-5-2`",
       "`pre_nm_55-6-1`","`pre_nm_56-1-3`","`pre_nm_56-3-2`","`pre_nm_63-1-1`")

pre_nm_41_13_1 <- lapply(`pre_nm_41-13-1`,function(x){ mutate(x,FB=(max.Hz-min.Hz)/1000)})
post_nm_41_13_1 <- lapply(`post_nm_41-13-1`,function(x){ mutate(x,FB=(max.Hz-min.Hz)/1000)})
pre_nm_49_12_1 <- lapply(`pre_nm_49-12-1`,function(x){ mutate(x,FB=(max.Hz-min.Hz)/1000)})
post_nm_49_12_1 <- lapply(`post_nm_49-12-1`,function(x){ mutate(x,FB=(max.Hz-min.Hz)/1000)})
pre_nm_53_2_1 <- lapply(`pre_nm_53-2-1`,function(x){ mutate(x,FB=(max.Hz-min.Hz)/1000)})
post_nm_53_2_1 <- lapply(`post_nm_53-2-1`,function(x){ mutate(x,FB=(max.Hz-min.Hz)/1000)})
pre_nm_54_2_2 <- lapply(`pre_nm_54-2-2`,function(x){ mutate(x,FB=(max.Hz-min.Hz)/1000)})
post_nm_54_2_2 <- lapply(`post_nm_54-2-2`,function(x){ mutate(x,FB=(max.Hz-min.Hz)/1000)})
pre_nm_55_4_3 <- lapply(`pre_nm_55-4-3`,function(x){ mutate(x,FB=(max.Hz-min.Hz)/1000)})
post_nm_55_4_3 <- lapply(`post_nm_55-4-3`,function(x){ mutate(x,FB=(max.Hz-min.Hz)/1000)})
pre_nm_55_5_2 <- lapply(`pre_nm_55-5-2`,function(x){ mutate(x,FB=(max.Hz-min.Hz)/1000)})
post_nm_55_5_2 <- lapply(`post_nm_55-5-2`,function(x){ mutate(x,FB=(max.Hz-min.Hz)/1000)})
pre_nm_55_6_1 <- lapply(`pre_nm_55-6-1`,function(x){ mutate(x,FB=(max.Hz-min.Hz)/1000)})
post_nm_55_6_1 <- lapply(`post_nm_55-6-1`,function(x){ mutate(x,FB=(max.Hz-min.Hz)/1000)})
pre_nm_56_1_3 <- lapply(`pre_nm_56-1-3`,function(x){ mutate(x,FB=(max.Hz-min.Hz)/1000)})
post_nm_56_1_3 <- lapply(`post_nm_56-1-3`,function(x){ mutate(x,FB=(max.Hz-min.Hz)/1000)})
pre_nm_56_3_2 <- lapply(`pre_nm_56-3-2`,function(x){ mutate(x,FB=(max.Hz-min.Hz)/1000)})
post_nm_56_3_2 <- lapply(`post_nm_56-3-2`,function(x){ mutate(x,FB=(max.Hz-min.Hz)/1000)})
pre_nm_63_1_1 <- lapply(`pre_nm_63-1-1`,function(x){ mutate(x,FB=(max.Hz-min.Hz)/1000)})
post_nm_63_1_1 <- lapply(`post_nm_63-1-1`,function(x){ mutate(x,FB=(max.Hz-min.Hz)/1000)})

# add in timepoint as a variable
pre_nm_41_13_1 <- lapply(pre_nm_41_13_1,function(x){ mutate(x,Timepoint="pre")})
post_nm_41_13_1 <- lapply(post_nm_41_13_1,function(x){ mutate(x,Timepoint="post")})
pre_nm_49_12_1 <- lapply(pre_nm_49_12_1,function(x){ mutate(x,Timepoint="pre")})
post_nm_49_12_1 <- lapply(post_nm_49_12_1,function(x){ mutate(x,Timepoint="post")})
pre_nm_53_2_1 <- lapply(pre_nm_53_2_1,function(x){ mutate(x,Timepoint="pre")})
post_nm_53_2_1 <- lapply(post_nm_53_2_1,function(x){ mutate(x,Timepoint="post")})
pre_nm_54_2_2 <- lapply(pre_nm_54_2_2,function(x){ mutate(x,Timepoint="pre")})
post_nm_54_2_2 <- lapply(post_nm_54_2_2,function(x){ mutate(x,Timepoint="post")})
pre_nm_55_4_3 <- lapply(pre_nm_55_4_3,function(x){ mutate(x,Timepoint="pre")})
post_nm_55_4_3 <- lapply(post_nm_55_4_3,function(x){ mutate(x,Timepoint="post")})
pre_nm_55_5_2 <- lapply(pre_nm_55_5_2,function(x){ mutate(x,Timepoint="pre")})
post_nm_55_5_2 <- lapply(post_nm_55_5_2,function(x){ mutate(x,Timepoint="post")})
pre_nm_55_6_1 <- lapply(pre_nm_55_6_1,function(x){ mutate(x,Timepoint="pre")})
post_nm_55_6_1 <- lapply(post_nm_55_6_1,function(x){ mutate(x,Timepoint="post")})
pre_nm_56_1_3 <- lapply(pre_nm_56_1_3,function(x){ mutate(x,Timepoint="pre")})
post_nm_56_1_3 <- lapply(post_nm_56_1_3,function(x){ mutate(x,Timepoint="post")})
pre_nm_56_3_2 <- lapply(pre_nm_56_3_2,function(x){ mutate(x,Timepoint="pre")})
post_nm_56_3_2 <- lapply(post_nm_56_3_2,function(x){ mutate(x,Timepoint="post")})
pre_nm_63_1_1 <- lapply(pre_nm_63_1_1,function(x){ mutate(x,Timepoint="pre")})
post_nm_63_1_1 <- lapply(post_nm_63_1_1,function(x){ mutate(x,Timepoint="post")})

# Combine pre- and post-surgery tables for each individual
nm_41_13_1 <- c(pre_nm_41_13_1,post_nm_41_13_1)
nm_49_12_1 <- c(pre_nm_49_12_1,post_nm_49_12_1)
nm_53_2_1 <- c(pre_nm_53_2_1,post_nm_53_2_1)
nm_54_2_2 <- c(pre_nm_54_2_2,post_nm_54_2_2)
nm_55_4_3 <- c(pre_nm_55_4_3,post_nm_55_4_3)
nm_55_5_2 <- c(pre_nm_55_5_2,post_nm_55_5_2)
nm_55_6_1 <- c(pre_nm_55_6_1,post_nm_55_6_1)
nm_56_1_3 <- c(pre_nm_56_1_3,post_nm_56_1_3)
nm_56_3_2 <- c(pre_nm_56_3_2,post_nm_56_3_2)
nm_63_1_1 <- c(pre_nm_63_1_1,post_nm_63_1_1)

# calculate arc length:chord length and add as a variable
nm_41_13_1 <- lapply(nm_41_13_1,function(x){ mutate(x,ACratio=arcL/chordL)})
nm_49_12_1 <- lapply(nm_49_12_1,function(x){ mutate(x,ACratio=arcL/chordL)})
nm_53_2_1 <- lapply(nm_53_2_1,function(x){ mutate(x,ACratio=arcL/chordL)})
nm_54_2_2 <- lapply(nm_54_2_2,function(x){ mutate(x,ACratio=arcL/chordL)})
nm_55_4_3 <- lapply(nm_55_4_3,function(x){ mutate(x,ACratio=arcL/chordL)})
nm_55_5_2 <- lapply(nm_55_5_2,function(x){ mutate(x,ACratio=arcL/chordL)})
nm_55_6_1 <- lapply(nm_55_6_1,function(x){ mutate(x,ACratio=arcL/chordL)})
nm_56_1_3 <- lapply(nm_56_1_3,function(x){ mutate(x,ACratio=arcL/chordL)})
nm_56_3_2 <- lapply(nm_56_3_2,function(x){ mutate(x,ACratio=arcL/chordL)})
nm_63_1_1 <- lapply(nm_63_1_1,function(x){ mutate(x,ACratio=arcL/chordL)})

# calculate theoretical max arc length possible to be able to detect notes where arc length measurement didn't work well
th.max_AL <- c()
th.max_ACr <- c()
for (i in 1:length(nm_63_1_1)){
  for (j in 1:nrow(nm_63_1_1[[i]])){
    a <- nm_63_1_1[[i]]$max.Hz[j] - nm_63_1_1[[i]]$min.Hz[j]
    b <- nm_63_1_1[[i]]$note.durs[j]/1000
    c <- sqrt((a^2)+(b^2))
    mAL <- a + b
    mACr <- (a+b)/c
    th.max_AL <- append(th.max_AL,mAL)
    th.max_ACr <- append(th.max_ACr,mACr)
  }
  nm_63_1_1[[i]]$th.maxAL <- th.max_AL
  nm_63_1_1[[i]]$th.maxACr <- th.max_ACr
  th.max_AL <- c()
  th.max_ACr <- c()
}

# visualize how many notes are above and below theoretical max arc length value
g1 <- ggplot()
for (i in 1:length(nm_41_13_1)) {
  g1 <- g1 + geom_point(data = nm_41_13_1[[i]], aes(x=th.maxAL,y=arcL),alpha=0.5,col="gray")
}
g1 + geom_abline(intercept = 0, slope = 1) +
  xlab("theoretical max arc length") + ylab("calculated arc length") +
  ggtitle("41-13-1") + theme_classic()

g2 <- ggplot()
for (i in 1:length(nm_49_12_1)) {
  g2 <- g2 + geom_point(data = nm_49_12_1[[i]], aes(x=th.maxAL,y=arcL),alpha=0.5,col="gray")
}
g2 + geom_abline(intercept = 0, slope = 1) +
  xlab("theoretical max arc length") + ylab("calculated arc length") +
  ggtitle("49-12-1") + theme_classic()

g3 <- ggplot()
for (i in 1:length(nm_53_2_1)) {
  g3 <- g3 + geom_point(data = nm_53_2_1[[i]], aes(x=th.maxAL,y=arcL),alpha=0.5,col="gray")
}
g3 + geom_abline(intercept = 0, slope = 1) +
  xlab("theoretical max arc length") + ylab("calculated arc length") +
  ggtitle("53-2-1") + theme_classic()

g4 <- ggplot()
for (i in 1:length(nm_54_2_2)) {
  g4 <- g4 + geom_point(data = nm_54_2_2[[i]], aes(x=th.maxAL,y=arcL),alpha=0.5,col="gray")
}
g4 + geom_abline(intercept = 0, slope = 1) +
  xlab("theoretical max arc length") + ylab("calculated arc length") +
  ggtitle("54-2-2") + theme_classic()

g5 <- ggplot()
for (i in 1:length(nm_55_4_3)) {
  g5 <- g5 + geom_point(data = nm_55_4_3[[i]], aes(x=th.maxAL,y=arcL),alpha=0.5,col="gray")
}
g5 + geom_abline(intercept = 0, slope = 1) +
  xlab("theoretical max arc length") + ylab("calculated arc length") +
  ggtitle("55-4-3") + theme_classic()

g6 <- ggplot()
for (i in 1:length(nm_55_5_2)) {
  g6 <- g6 + geom_point(data = nm_55_5_2[[i]], aes(x=th.maxAL,y=arcL),alpha=0.5,col="gray")
}
g6 + geom_abline(intercept = 0, slope = 1) +
  xlab("theoretical max arc length") + ylab("calculated arc length") +
  ggtitle("55-5-2") + theme_classic()

g7 <- ggplot()
for (i in 1:length(nm_55_6_1)) {
  g7 <- g7 + geom_point(data = nm_55_6_1[[i]], aes(x=th.maxAL,y=arcL),alpha=0.5,col="gray")
}
g7 + geom_abline(intercept = 0, slope = 1) +
  xlab("theoretical max arc length") + ylab("calculated arc length") +
  ggtitle("55-6-1") + theme_classic()

g8 <- ggplot()
for (i in 1:length(nm_56_1_3)) {
  g8 <- g8 + geom_point(data = nm_56_1_3[[i]], aes(x=th.maxAL,y=arcL),alpha=0.5,col="gray")
}
g8 + geom_abline(intercept = 0, slope = 1) +
  xlab("theoretical max arc length") + ylab("calculated arc length") +
  ggtitle("56-1-3") + theme_classic()

g9 <- ggplot()
for (i in 1:length(nm_56_3_2)) {
  g9 <- g9 + geom_point(data = nm_56_3_2[[i]], aes(x=th.maxAL,y=arcL),alpha=0.5,col="gray")
}
g9 + geom_abline(intercept = 0, slope = 1) +
  xlab("theoretical max arc length") + ylab("calculated arc length") +
  ggtitle("56-3-2") + theme_classic()

g10 <- ggplot()
for (i in 1:length(nm_63_1_1)) {
  g10 <- g10 + geom_point(data = nm_63_1_1[[i]], aes(x=th.maxAL,y=arcL),alpha=0.5,col="gray")
}
g10 + geom_abline(intercept = 0, slope = 1) +
  xlab("theoretical max arc length") + ylab("calculated arc length") +
  ggtitle("63-1-1") + theme_classic()

# only keep notes that are ACr < = theoretical max
reduced_nm_41_13_1 <- nm_41_13_1 # new copy of list
reduced_nm_49_12_1 <- nm_49_12_1
reduced_nm_53_2_1 <- nm_53_2_1
reduced_nm_54_2_2 <- nm_54_2_2
reduced_nm_55_4_3 <- nm_55_4_3
reduced_nm_55_5_2 <- nm_55_5_2
reduced_nm_55_6_1 <- nm_55_6_1
reduced_nm_56_1_3 <- nm_56_1_3
reduced_nm_56_3_2 <- nm_56_3_2
reduced_nm_63_1_1 <- nm_63_1_1

bad <- c()
for (i in 1:length(reduced_nm_63_1_1)){
  if (nrow(reduced_nm_63_1_1[[i]]) > 0){
    for (j in 1:nrow(reduced_nm_63_1_1[[i]])){
      if (reduced_nm_63_1_1[[i]]$ACratio[j] > reduced_nm_63_1_1[[i]]$th.maxACr[j]){
        bad <- append(bad,j)
      }
    }
    if (!is.null(bad)){
      reduced_nm_63_1_1[[i]] <- reduced_nm_63_1_1[[i]][-bad,]
      bad <- c()
    }
  }
}

# only keep notes that are = or < theoretical max AL
redAL_nm_41_13_1 <- nm_41_13_1 # new copy of list
redAL_nm_49_12_1 <- nm_49_12_1
redAL_nm_53_2_1 <- nm_53_2_1
redAL_nm_54_2_2 <- nm_54_2_2
redAL_nm_55_4_3 <- nm_55_4_3
redAL_nm_55_5_2 <- nm_55_5_2
redAL_nm_55_6_1 <- nm_55_6_1
redAL_nm_56_1_3 <- nm_56_1_3
redAL_nm_56_3_2 <- nm_56_3_2
redAL_nm_63_1_1 <- nm_63_1_1

bad <- c()
for (i in 1:length(redAL_nm_63_1_1)){
  if (nrow(redAL_nm_63_1_1[[i]]) > 0){
    for (j in 1:nrow(redAL_nm_63_1_1[[i]])){
      if (redAL_nm_63_1_1[[i]]$arcL[j] > redAL_nm_63_1_1[[i]]$th.maxAL[j]){
        bad <- append(bad,j)
      }
    }
    if (!is.null(bad)){
      redAL_nm_63_1_1[[i]] <- redAL_nm_63_1_1[[i]][-bad,]
      bad <- c()
    }
  }
}

save(redAL_nm_41_13_1,redAL_nm_49_12_1,redAL_nm_53_2_1,redAL_nm_54_2_2,redAL_nm_55_4_3,
     redAL_nm_55_5_2,redAL_nm_55_6_1,redAL_nm_56_1_3,redAL_nm_56_3_2,redAL_nm_63_1_1,
     file="reducedAL_lists.RData")

g1 <- ggplot()
for (i in 1:length(redAL_nm_41_13_1)) {
  g1 <- g1 + geom_point(data = redAL_nm_41_13_1[[i]], 
                        aes(x=th.maxAL,y=arcL),alpha=0.5) + 
    geom_abline(slope=1,intercept = 0)
}
g1 + xlab("theoretical max arc length") + ylab("arc length") + 
  ggtitle("41-13-1 reduced dataset") + theme_classic()

lengths <- data.frame(ID = c("41-13-1","41-13-1","49-12-1","49-12-1","53-2-1","53-2-1","54-2-2","54-2-2",
                             "55-4-3","55-4-3","55-5-2","55-5-2","55-6-1","55-6-1","56-1-3","56-1-3",
                             "56-3-2","56-3-2","63-1-1","63-1-1"), 
                      Timepoint = rep(c("pre","post"),10),
                      length = c(length(pre_nm_41_13_1), length(post_nm_41_13_1),
                                 length(pre_nm_49_12_1), length(post_nm_49_12_1),
                                 length(pre_nm_53_2_1), length(post_nm_53_2_1),
                                 length(pre_nm_54_2_2), length(post_nm_54_2_2),
                                 length(pre_nm_55_4_3), length(post_nm_55_4_3),
                                 length(pre_nm_55_5_2), length(post_nm_55_5_2),
                                 length(pre_nm_55_6_1), length(post_nm_55_6_1),
                                 length(pre_nm_56_1_3), length(post_nm_56_1_3),
                                 length(pre_nm_56_3_2), length(post_nm_56_3_2),
                                 length(pre_nm_63_1_1), length(post_nm_63_1_1)))

save(nm_41_13_1,nm_49_12_1,nm_53_2_1,nm_54_2_2,nm_55_4_3,
     nm_55_5_2,nm_55_6_1,nm_56_1_3,nm_56_3_2,nm_63_1_1,
     redAL_nm_41_13_1,redAL_nm_49_12_1,redAL_nm_53_2_1,
     redAL_nm_54_2_2,redAL_nm_55_4_3,redAL_nm_55_5_2,
     redAL_nm_55_6_1,redAL_nm_56_1_3,redAL_nm_56_3_2,
     redAL_nm_63_1_1,lengths,file="new_note_matrices.RData")
