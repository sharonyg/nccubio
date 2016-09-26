pam1<-read.table("D:\NCCU\Bioinformations\hw1\pam.text", sep=",")
colnames(pam1)<-c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
rownames(pam1)<-c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")

data2 <- as.matrix(pam1)
data3 <- data2/10000

mat <- round(mpower(data3, 250), digits=3)
pam250 <- round(mat*100, digits = 0)
pam250

write.table(pam250, file = "pam250.txt", sep = "\t", row.names = T, col.names = T)
#write.matrix(pam250,'pam250.txt',sep = "")
