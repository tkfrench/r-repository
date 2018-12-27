
my2by2.test <- function(num1,den1,num2,den2,muex=T){
  if(muex==FALSE){
      t<-matrix(c(num1, den1-num1, num2-num1, den2-num2-(den1-num1)),nrow=2)
  } else{
    t<-matrix(c(num1, den1-num1, num2, den2-num2), nrow=2)
  }
  print( (num1+num2)*((den1-num1)+(den2-num2)) )
  # test if too big for Fisher test
  if(((num1+num2)*((den1-num1)+(den2-num2)))< 2^31 - 1)
    {fisher.test(t)}
  else{chisq.test(t)}
}

# Given Typical Proportion counts from SQL as:
#   Hosp: (a) num=5    (b) den=100
#   sys:  (c) num=100, (d) den=1000
#  
#    test:  hosp  v system      Mutually exl. =TRUE
#         ------------------+   +-----------+
#      w  |   5   |   100   |   | a   | c   | 
#         +-----------------+   +-----------+
#     W/o |   95  |   900   |   | b-a | d-c |  
#         +-----------------+   +-----------+
#     tot |   100 |   1000  |   | b   | d   |
my2by2.test(5,100,100,1000,muex=T)
fisher.test(matrix(c(5,95,100,900),nrow=2)) # should be the same

#    test:   hosp  v non-hosp   Mutually exl. =FALSE
#         ------------------+   +-----------------+
#      w  |   5   |    95   |   | a   | c-a       | 
#         +-----------------+   +-----------------+
#     W/o |   95  |   805   |   | b-a | d-c-(b-a) |  
#         +-----------------+   +-----+-----------+
#     tot |   100 | 1000-100    | b   | d- b  
my2by2.test(5,100,100,1000,muex=F)
fisher.test(matrix(c(5,95,95,805),nrow=2)) # should be the same

my2by2.test(240000,500000,250000,500000,muex=T)
fisher.test(matrix(c(250000,500000,250000,500000),nrow=2)) 
