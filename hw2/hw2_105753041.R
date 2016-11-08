library(Biostrings)
library(seqinr)
library(rnaseqWrapper)

#讀進序列
input_seq <- read.alignment("test.fasta", format="fasta",  forceToLower = FALSE)
input_seq$seq
#選擇需要的matrix(PAM100 or PAM250...etc.)
inputMatrix <- read.table(file.choose())

#創建新的dataframe
df <- data.frame(matrix(ncol = 58, nrow = 61))
names(df)[1:58] <- sprintf("%d", 1:58)

#分割兩條sequence
col_string <- strsplit(input_seq$seq[1], "")
col_string
row_string <- strsplit(input_seq$seq[2], "")
row_string

#list中加上gap
col_string_matrix <- append(col_string[[1]], "X.")
row_string_matrix <- append(row_string[[1]], "*")
#把gap搬到最前面
a <- length(col_string_matrix)
b <- length(row_string_matrix)
while(a>1)
{
  col_string_matrix[a] <- col_string_matrix[a-1]
  a <- a-1
}
while(b>1)
{
  row_string_matrix[b] <- row_string_matrix[b-1]
  b <- b-1
}
col_string_matrix[1] <- "*"
row_string_matrix[1] <- "*"
col_string_matrix
row_string_matrix

#填滿第一行第一列(gap)
df[1,1] <- 0
gap_open <- -10
for (x in c(1:(ncol(df)-1)))
{
  df[1,x+1] <- df[1,x] + gap_open
  x <- x+1
}
for (y in c(1:(nrow(df)-1)))
{
  df[y+1,1] <- df[y,1] + gap_open
  y <- y+1
}

#score matrix
for(i in c(2:(nrow(df))))
{
  for(j in c(2:(ncol(df))))
  {
    up_score <- df[i-1,j] + gap_open
    left_score <- df[i,j-1] + gap_open
    diag_score <- df[i-1,j-1] + inputMatrix[row_string_matrix[j], col_string_matrix[i]]
    max_score <- max(up_score, left_score, diag_score)
    df[i,j] <- max_score
    j <- j+1
  }
  i <- i+1
}

#找最大值
maxIndex <- which(df == max(df), arr.ind = TRUE)
max_index <- apply(maxIndex,2,max)
i <- max_index[1]
j <- max_index[2]
df[i,j]

#trace back & alignment
seq_aln <- list(seq1=c(), seq2=c())
seq_aln

while(i>1 && j>1)
{
  trace_max <- max(df[i,j-1], df[i-1,j-1], df[i-1,j])
  if(trace_max == df[i-1,j-1])#走斜角(Substitution)
  {
    seq_aln$seq1 <- append(seq_aln$seq1, col_string_matrix[j])
    seq_aln$seq2 <- append(seq_aln$seq2, row_string_matrix[i])
    i <- i-1
    j <- j-1
  }
  else if(trace_max == df[i-1,j])#走上面(Insertion)
  {
    seq_aln$seq1 <- append(seq_aln$seq1, "-")
    seq_aln$seq2 <- append(seq_aln$seq2, row_string_matrix[i])
    i <- i-1
  }
  else if(trace_max == df[i,j-1])#走左邊(Deletion)
  {
    seq_aln$seq1 <- append(seq_aln$seq1, col_string_matrix[j])
    seq_aln$seq2 <- append(seq_aln$seq2, "-")
    j <- j-1
  }
  else
    break
}

seq_aln$seq1 <- paste("", rev(seq_aln$seq1), sep="", collapse="")
seq_aln$seq2 <- paste("", rev(seq_aln$seq2), sep="", collapse="")

#output
write.fasta(sequences=seq_aln, names=names(seq_aln),file.out="result.fasta")

