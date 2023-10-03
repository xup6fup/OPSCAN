
library(magrittr)
library(rtf)

# 0. Settings

intrested_var <- c('ISCD indication', 'total_CXR_group', 'positive_prop_group',
                   'GENDER', 'AGE_group', 'BMI_group',
                   'Fracture', 'SecondaryOsteoporosis', 
                   'RheumatoidArthritis', 'Steroid')

data_path <-  './data/RCT analysis/RCT data.RData'
table_path <- './result/Table 02.doc'

# 0. Functions

.ROUND <- function(X,digits=0) {
  x=rep(NA,length(X))
  for (a in 1:length(X)) {
    if (is.na(X[a])) {x[a]=NA} else {
      x[a]=formatC(X[a],format="f",digits=digits)
    }
  }
  return(x)
}

.ROUND.p <- function(X,p.digits=3) {
  p=rep(NA,length(X))
  for (a in 1:length(X)) {
    if (!is.na(X[a])) {
      if (p.digits==0) {p[a]="<1"} else {
        p[a]=as.character(round(X[a],p.digits))
        if (p[a]=="0") {
          if (p.digits==1) {p[a]="<0.1"} else {
            p[a]="<0."
            for (l in 1:(p.digits-1)) {
              p[a]=paste0(p[a],0)
            }
            p[a]=paste0(p[a],1)  
          }
        } else {
          p[a]=.ROUND(as.numeric(X[a]),p.digits)
        }
      }
    }
  }
  return(p)
}

Table1 <- function(X, Y.matrix,digits=1,digits.per=1,p.digits=3,Factor=NULL,x.name="Group",SD=TRUE,nonparametric=TRUE,Transpose=FALSE) {
  if (is.null(Factor)) {
    Factor=numeric()
    for (k in 1:ncol(Y.matrix)) {
      if (is.factor(Y.matrix[,k])) {
        Factor=c(Factor,1)
      } else {
        Factor=c(Factor,0)
      }   
    }
  }
  n.row=numeric()
  row.names=NULL
  for (k in 1:ncol(Y.matrix)) {
    if (Factor[k]==1) {
      n.row=c(n.row,1+length(levels(Y.matrix[,k])))
      row.names=c(row.names,colnames(Y.matrix)[k],paste0(colnames(Y.matrix)[k],":",levels(Y.matrix[,k])))
    } else {
      n.row=c(n.row,1)
      row.names=c(row.names,colnames(Y.matrix)[k])
    }   
  }
  Table=matrix("",nrow=sum(n.row),ncol=length(levels(factor(X)))+1)
  colnames(Table)=c(paste0(x.name,":",levels(factor(X))),"p-value")
  rownames(Table)=row.names
  position=cumsum(n.row)
  #######################################
  for (k in 1:ncol(Y.matrix)) {
    if (Factor[k]==0) {
      n.sample=numeric()     
      for (i in 1:length(levels(factor(X)))) {
        n.sample=c(n.sample,sum(is.na(Y.matrix[factor(X)==levels(factor(X))[i],k])==F))
        m=mean(Y.matrix[factor(X)==levels(factor(X))[i],k],na.rm=T)
        if (SD==TRUE) {
          s=sd(Y.matrix[factor(X)==levels(factor(X))[i],k],na.rm=T)
        } else {
          s=sd(Y.matrix[factor(X)==levels(factor(X))[i],k],na.rm=T)/sqrt(sum(!is.na(Y.matrix[factor(X)==levels(factor(X))[i],k])))
        }
        m=.ROUND(m,digits)
        s=.ROUND(s,digits)
        Table[position[k],i]=paste0(m,"±",s)
      }
      if (length(levels(factor(X)))>1) {
        p <- ''
        if (min(n.sample) > 0) {
          if (nonparametric==TRUE) {
            if (min(n.sample)>=25) {
              p=.ROUND.p(anova(lm(Y.matrix[,k]~factor(X)))$"Pr(>F)"[1],p.digits)
            } else {
              if (length(levels(factor(X)))>=3) {
                p=.ROUND.p(kruskal.test(Y.matrix[,k],factor(X))$p.value,p.digits)
              } else {
                p=.ROUND.p(wilcox.test(Y.matrix[,k]~factor(X),correct=FALSE)$p.value,p.digits)
              }
            }
            if (min(n.sample)<25) {p=paste0(p,"#")} 
          } else {
            p=.ROUND.p(anova(lm(Y.matrix[,k]~factor(X)))$"Pr(>F)"[1],p.digits)
          }
        }
        Table[position[k],length(levels(factor(X)))+1] = p
      }
    } else {
      for (i in 1:length(levels(factor(X)))) {
        for (j in 1:length(levels(Y.matrix[,k]))) {
          n=sum(Y.matrix[factor(X)==levels(factor(X))[i],k]==levels(Y.matrix[,k])[j],na.rm=T)
          prop=n/length(Y.matrix[factor(X)==levels(factor(X))[i]&is.na(Y.matrix[,k])==F&is.na(X)==F,k])
          prop=.ROUND(prop*100,digits.per)
          prop=paste0("(",prop,"%)")
          Table[position[k]-length(levels(Y.matrix[,k]))+j,i]=paste0(n,prop)
        }
      }
      if (length(levels(factor(X)))>1) {
        if (nonparametric==TRUE) {
          logic=tryCatch(chisq.test(table(Y.matrix[,k],factor(X))),error=function(e) e, warning=function(w) w)
          sig.test=try(fisher.test(table(Y.matrix[,k],factor(X))),silent=T)
          if (is(logic,"warning")&is(sig.test)!="try-error") {
            p=.ROUND.p(sig.test$p.value,p.digits)
          } else {
            p=.ROUND.p(chisq.test(table(Y.matrix[,k],factor(X)),correct=FALSE)$p.value,p.digits)
          }
          if (is(logic,"warning")&is(sig.test)!="try-error") {p=paste0(p,"#")}
        } else {
          p=.ROUND.p(chisq.test(table(Y.matrix[,k],factor(X)),correct=FALSE)$p.value,p.digits)
        }
        Table[position[k]-length(levels(Y.matrix[,k])),length(levels(factor(X)))+1]=p
      }
    }
  }
  Table.Final=cbind(rownames(Table),Table)
  colnames(Table.Final)[1]="Variable"
  rownames(Table.Final)=1:nrow(Table)
  if (length(levels(factor(X)))>1) {return(Table.Final)} else {
    Table.Final[,3]=""
    colnames(Table.Final)[2]="Total"
    return(Table.Final[,-3])
  }
}

Table2doc <- function (table_list, filename = 'test.doc') {
  
  rtffile <- RTF(filename, width = 8.5, height = 11,font.size = 10, omi = c(0.5, 0.5, 0.5, 0.5))
  
  for (j in 1:length(table_list)) {
    
    addHeader(rtffile, paste0("Table ", j))
    ALIGN <- c("L", rep("C", ncol(table_list[[j]]) - 1))
    new_table <- table_list[[j]]
    addTable(rtffile, new_table, NA.string = "", row.names = FALSE, col.justify = ALIGN, header.col.justify = ALIGN)
    addNewLine(rtffile)
    
  }
  
  done(rtffile)
  
}

# 1. Load data

load(data_path)

follow_data <- follow_data[follow_data[,'Group'] %in% c('Screening', 'Control'),]
follow_data[,'Group'] <- as.character(follow_data[,'Group']) %>% factor(., levels = c('Screening', 'Control'))

follow_data[,'ISCD indication'] <- factor((follow_data[,'GENDER'] %in% 'male' & follow_data[,'AGE'] >= 70) | (follow_data[,'GENDER'] %in% 'female' & follow_data[,'AGE'] >= 65))
levels(follow_data[,'ISCD indication']) <- c('Not meeting ISCD indication', 'Meeting ISCD indication')

# 2. Data

data_list <- list()

data_list[[1]] <- follow_data[follow_data[,'Group'] %in% 'Screening',]
data_list[[1]][,'y'] <- (data_list[[1]][,'DXA_type'] %in% 2) + 0L

data_list[[2]] <- follow_data[follow_data[,'Group'] %in% 'Screening' & follow_data[,'DXA_type'] %in% 2,]
data_list[[2]][,'y'] <- (data_list[[2]][,'OP_type'] %in% 2) + 0L

data_list[[3]] <- follow_data[follow_data[,'Group'] %in% 'Screening' & follow_data[,'DXA_type'] %in% 2 & follow_data[,'OP_type'] %in% 2,]
data_list[[3]][,'y'] <- (data_list[[3]][,'AOM_type'] %in% 1) + 0L

# 3. Build table

table_list <- list()

for (i in 1:length(data_list)) {
  
  sub_table <- NULL
  
  for (j in 1:length(intrested_var)) {
    
    count_tab <- table(data_list[[i]][,c(intrested_var[j], 'y')])
    prop_tab <- prop.table(count_tab, 1)
    
    current_table <- matrix('', nrow = nrow(count_tab) + 1, ncol = 3)
    
    current_table[1,1] <- intrested_var[j]
    current_table[-1,1] <- rownames(count_tab)
    current_table[-1,2] <- paste0(count_tab[,'1'], '/', apply(count_tab, 1, sum), ' (', formatC(prop_tab[,'1'] * 100, 1, format = 'f'), '%)')
    
    chisq_result <- try(chisq.test(count_tab), silent = TRUE)
    
    if (!'try-error' %in% class(chisq_result)) {
      
      p_val <- paste0('p = ', formatC(chisq_result[['p.value']], 3, format = 'f'))
      if (p_val %in% 'p = 0.000') {p_val <- 'p < 0.001'}
      current_table[1,3] <- p_val
      
    }
    
    sub_table <- rbind(sub_table, current_table)
    
  }
  
  if (!i %in% 1) {sub_table <- sub_table[,-1,drop=FALSE]}
  
  table_list[[i]] <- sub_table

}

# 4. Write out

Table2doc(table_list = list(do.call('cbind', table_list)), filename = table_path)
