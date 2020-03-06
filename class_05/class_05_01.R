#' ---
#' title: "Class 05: Data Visualization and grahpics in R"
#' author: "Brandon Jones"
#' date: "Jan 24th, 2020"
#' ---
#' 
#Class 05
#Data Visualization and grahpics in R
plot(1:5, col="blue", typ="o")
weight_chart <- read.table("bimm143_05_rstats/weight_chart.txt",header=TRUE)
Mouse <- read.table("bimm143_05_rstats/feature_counts.txt",
                    header=TRUE, sep="\t")

plot(weight_chart$Age, weight_chart$Weight, type="p", pch=15)
plot(weight_chart$Age, weight_chart$Weight, type="p", pch=15, cex=1.5)
plot(weight_chart$Age, weight_chart$Weight, type="o", pch=15, col="blue", 
     cex=1.2, main = "Baby weight with age",
     xlab="Age(months)",ylab="Weight(kg)",xlim=c(0,9),ylim=c(2,10))
Mouse <- read.table("bimm143_05_rstats/feature_counts.txt",header=TRUE, sep="\t")

gender <- read.delim("bimm143_05_rstats/male_female_counts.txt", header=TRUE)
barplot(gender$Count, names.arg = gender$Sample, col=rainbow (nrow(gender)) )


par(las=1, mar=c(3,11,2,1))
barplot(Mouse$Count,horiz=TRUE, names.arg = Mouse$Feature, ylab="Feature",
        main="Mouse Genome Feature Count")





