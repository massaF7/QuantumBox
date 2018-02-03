date
    "x_1" = c(1,0,1,0,1,0,1,0,1,0,1,0,1,0) # 入力
    "x_2" = c(0,1,1,0,0,0,0,1,1,1,0,0,0,1)
    "x_3" = c(1,0,1,1,1,1,0,1,0,1,0,1,0,1)
    "x_4" = c(0,1,1,1,1,0,0,0,0,0,1,0,0,1)
    "x_5" = c(1,0,0,0,0,0,1,1,1,1,0,1,0,0)
    "x_6" = c(0,1,0,1,0,1,0,1,1,0,1,1,0,0)
    "x_7" = c(1,0,0,1,1,0,1,1,0,1,0,1,1,0)

    "w_1" = w # 入力
    "w_2" = w
    "w_3" = w
    "w_4" = w
    "w_5" = w
    "w_6" = w
    "w_7" = w

    "AoCw_1" = w # 入力
    "AoCw_2" = w
    "AoCw_3" = w
    "AoCw_4" = w
    "AoCw_5" = w
    "AoCw_6" = w
    "AoCw_7" = w

t = c(0,1,0,1,0,0,1,0,0,1,0,0,0,1)# 教師信号
w = 12*(2*2-3)/50 # シナプス結合荷重
h = (12/50)+0.1# 閾値
η = c(0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0)

# シナプス結合荷重、閾値の初期値を設定
# 入力パターンの中から1つを選び、それに対応する出力を計算
# 入力パターンに対応する教師信号との差を計算し、シナプス結合荷重、閾値を修正
# 誤りが０となれば終了、それ以外は繰り返し
while η=η+1 η=8{
    while t=t+1 t=14{
        f(u)
            if Z>0,z<-1,z<-0{
                if z-t=0, End,{
                    w_1 = w + x_1 * η * (t-z)
                    w_2 = w + x_2 * η * (t-z)
                    w_3 = w + x_3 * η * (t-z)
                    w_4 = w + x_4 * η * (t-z)
                    w_5 = w + x_5 * η * (t-z)
                    w_6 = w + x_6 * η * (t-z)
                    w_7 = w + x_7 * η * (t-z)
                    h=-η(t-z)
                }
            }
    }
}

E_All<-function(m,l,s,t,N,B,value=TRUE){
    A=B/N
    Rt=c()                           #generate infinite dimensional vector Rq
    Rq=c()
    Rqb=c()
    Rq2=c()
    Rq3=c()
    Rq4=c()
    CRq=c()
    CRq2=c()
    CRq3=c()
    CRq4=c()
    FullCRq=c()
    CRqb=c()
    M=c()
    M2=c()
    M3=c()
    M4=c()
    Mb=c()
    Rt2=c()
    Rt3=c()
    Rt4=c()
    Rtb=c()
    FullRt=c()
    Rttest=c()
    Rttest2=c()

    for(k in 1:N){
        b=A*k

        Rq[k]=ask(m,b)           ##cutoff-rate by ipop
        Rqb[k]=BanOpt(m,b)    #cutoff-rate by Ban's formula
        Rq2[k]=ask(l,b)           ##cutoff-rate for binary states
    		Rq3[k]=ask(s,b)	##
    		Rq4[k]=ask(t,b)	##

        M[k]=askene(m,b)        ##energy for optimal parameter obtained by ipop
        Mb[k]=BanOptEne(m,b) #energy for optimal parameter obtained by Ban's formula
        M2[k]=askene(l,b)        ##energy for binary states
    		M3[k]=askene(s,b)	##
    		M4[k]=askene(t,b)	##

        CRq[k]=Mu(M[k])   ##continuous cutoff-rate for energy M[k]
        CRq2[k]=Mu(M2[k]) ##continuous cutoff-rate for energy M2[k]
    		CRq3[k]=Mu(M3[k]) ##
    		CRq4[k]=Mu(M4[k]) ##
        CRqb[k]=Mu(Mb[k]) #continuous cutoff-rate for energy Mb[k]
        FullCRq[k]=MuFull(M[k])

        Rt[k]=Rq[k]/CRq[k] ##
        Rt2[k]=Rq2[k]/CRq2[k] ##
    		Rt3[k]=Rq3[k]/CRq3[k] ##
    		Rt4[k]=Rq4[k]/CRq4[k] ##
        Rtb[k]=Rqb[k]/CRqb[k]
        FullRt[k]=CRq[k]/FullCRq[k]
        Rttest[k]=Rqb[k]/Rq[k]
        Rttest2[k]=Mb[k]/M[k]
        cat(sprintf("%f,%e:%e/%e=%e¥n",b,M[k],Rq[k],CRq[k],Rt[k]))

    }

    Ubx=M[N]
    Uby=CRq[N]
    min_x=0
    min_y=0
    max_x=200
    max_y=4
    if(value){
            plot(M,Rq,lwd=1,lty=2,xlim=c(min_x,max_x),ylim=c(min_y,max_y),xlab="",ylab="",type="l")　
                par(new=TRUE)
            plot(M2,Rq2,lwd=1,lty=3,xlim=c(min_x,max_x),ylim=c(min_y,max_y),xlab="",ylab="",type="l",col="black")
                par(new=TRUE)
            plot(M3,Rq3,lwd=1,lty=4,xlim=c(min_x,max_x),ylim=c(min_y,max_y),xlab="",ylab="",type="l",col="black")
                par(new=TRUE)
            plot(M4,Rq4,lwd=1,lty=1,xlim=c(min_x,max_x),ylim=c(min_y,max_y),xlab="",ylab="",type="l",col="black")
        }else{
            lwds=c(2,1)
            ltys=c(2,1)
            plot(M,Rt,lwd=2,lty=2,xlim=c(min_x,max_x),ylim=c(min_y,max_y),xlab="",ylab="",type="l",col="black") ##
                par(new=TRUE) ##
            plot(M,Rt2,lwd=1,lty=3,xlim=c(min_x,max_x),ylim=c(min_y,max_y),xlab="",ylab="",type="l",col="black") ##
                par(new=TRUE)
            plot(M,Rt3,lwd=1,lty=4,xlim=c(min_x,max_x),ylim=c(min_y,max_y),xlab="",ylab="",type="l",col="black") ##
                par(new=TRUE)
            plot(M,Rt4,lwd=1,lty=1,xlim=c(min_x,max_x),ylim=c(min_y,max_y),xlab="",ylab="",type="l",col="black") ##
      }
    return(list(Rq,CRq,Rt,Rttest))
}

GetGraph<-function(){
m=5
N=100
B=50
l=10
s=20
t=30
pdf(file="yymmdd_ann.pdf")
     E_All(m,l,s,t,N,B,TRUE)
     dev.off()

}
