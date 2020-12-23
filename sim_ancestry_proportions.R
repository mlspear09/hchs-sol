##########
##
##  This is a Moran model style simulator, which, by definition, allows for overlapping generations.
##  Features include population growth (exponential), migration (with adjustable levels of migration & ancestry
##  patterns), ancestry-based assortative mating, and ancestry-based variability in fecundity.
##
##  The goal is to simulate Amerindigenous ancestry patterns seen in Mexican Americans, reported in Spear et al (2020).
##
##
## KEY FUNCTIONS:
##
##     simAI() - takes all parameters and generates an output table with nr=#sims nc=#generations
##
##     plotAIsim() - takes all parameters, runs simAI(), and generates stand-alone plot
##
##
##
##########

require(RColorBrewer)
pcol=brewer.pal(9, "Set1")
pcol=c(pcol,"black")
tmp=pcol[1];
pcol[1] = pcol[7];
pcol[7] = pcol[5];
pcol[5] = tmp;

## DEFINE PARAMETERS ##

NSIMS = 2; #total number of simulations to run
NGENS = 2; #total number of generations to simulate

N0 = 100; #initial pop size
G = 0; #pop growth rate, exponential

AI0 = 0.42; #initial AI ancestry proportion

M0 = 0; #initial fraction population constituted from migrants
MG = 0; #growth rate of migrant fraction
MAI0 = 0; #migrant AI parameter. 0=>same as domestic pop, 1=>migrants have AI=1.
MAIG = 0; #migrant AI growth rate. 0=>constant, large=>increasing AI% in migrant pool

FAI = 0; #fecundity AI parameter. 0=>uniform, 1=>higher fecundity for higher AI: rbeta(1, 1, 1/(1+FAI))

PSI = 0; #Assortative mating param between [0,1]: weighting from ~uniform (PSI=0) to N(Parent 1, 10)

EPS = 1e-8; #small number for rounding to 0.

## function to perform simulations. See above for parameter definitions. Defaults to SNM
simAI <- function(NSIMS=100, NGENS=2, N0=1000, G=0, AI0=0.42, M0=0, MG=0, MAI0=0, MAIG=0,
                  FAI=0, PSI=0, PLOT=0){
    ##include a status bar
    pb <- txtProgressBar(min = 0, max = NSIMS, style = 3)

    AIout = matrix(NA, nr=NSIMS, nc=NGENS+1); #output trajectory matrix
    if(PSI < 0){
        PSI=0;
    }
    if(PSI > 1){
        PSI = 1;
    }
    
    if(PLOT == 1){ #plot updates to default device
        plot(1,1,type='n',xlim=c(0,NGENS),ylim=c(0.37,0.6), xlab="Gen", ylab="AI anc")
    }
    ## START SIMULATION ##
    for(sim in 1:NSIMS){
        ##cat("working on sim ",sim,"/",NSIMS,"\n",sep="")
        
        ## initialization
        N = N0; # N is always population size at beginning of each generation
        M = M0
        AI = rbeta(N, 4*AI0/(1-AI0), 4) #initial population
        MAI = MAI0;
        AIout[sim, 1] = mean(AI);
        ## simulate generations
        for(gen in 1:NGENS){
            ##cat("\tworking on gen ",gen,"/",NGENS,"\n",sep="");
            NtF = round(N*(1+G)) # population size at end of this generation (final)
            MAIF = min(MAI*(1+MAIG), 1);
            MF = min(M*(1+MG), 1);
            
            ## each generation has NtF steps to reconstitute new population
            for(it in 1:NtF){
                Nt = round(N*(1+G)^(it/NtF)); # population size for this iteration
                maxP = ifelse(it<=N, N, it-2); #maximum parent ID to be defined
                MAIt = min(MAI*(1+MAIG)^(it/NtF), 1); # mean AI of migrant pop in this iteration
                Mt = M*(1+MG)^(it/NtF); # fraction of pop that is migrant in this iteration.
                ##cat("\tworking on it ",it,"/",NtF,"; Nt=",Nt,"\n",sep="");
                
                ## need to choose an individual to replace.
                ##  If there is exp growth, then we need to add new individuals.
                ##  So: for the first N steps (N = pop size at beginning of generation)
                ##  randomly sample an individual. Then sequentially add new individuals.

                if(it <= maxP){ #individual already exists, randomly sample one of them
                    i = sample(maxP,1) #randomly sample individual to replace
                }
                else{ #need to add new individual to population
                    i = maxP;
                }
                
                ## each individual is either a migrant (sample new individual) or not (randomly
                ## sample two parents).
                if(runif(1) <= Mt){ # migrant!
                    if(MAIt < EPS){
                        ## no difference between migrant and domestic AI%,
                        ##so random sample from domestic pop, randomly 1:(Nt-1)
                        AI[i] = AI[sample(maxP-1,1)]; 
                    }
                    else if(abs(MAIt-1) < EPS){ #complete AI ancestry (AI=1), so just set to 1
                        AI[i] = 1;
                    }
                    else{ #sample from shifted AI distribution

#########
### NEED TO FIX THIS TO SPECIFY ARBITRARY MAI
#########
                        t = MAIt + (1-MAIt)*mean(AI); #weighted average of current AI and 1.
                        AI[i] = rbeta(1, 4*t/(1-t), 4)
                    }
                }
                else{ # choose domestic parents
                    ## parent 1 chosen uniformly (if FAI==1) or biased toward higher AI (FAI>1)
                    if(FAI < EPS){ #no fecundity differences, random sample
                        P1 = sample(maxP-2,1); #randomly choose parent between 1:N
                    }
                    else{ #fecundity differences
                        ## draw P1 from (potentially) biased sample of domestic pop
                        ## since AI is not sorted, order(AI) gives AI rank of each individual, with
                        ## higher rank indicating higher AI%.
                        P1 = order(AI[1:(maxP-1)])[round((maxP-2)*rbeta(1,1,1/(1+FAI)))+1]; 
                    }
                    
                    ##choose parent 2.
                    ##  assortative mating, choose P2 index based on N(P1, sigma^2).
                    ##   sigma^2 is large (Nt) when PSI=0, and small (10) when PSI=1.
                    ##   sigma^2 is linear between these endpoints
                    ##   Note that P2 is truncated Normal so it is between [1,Nt-1].
                    ##   sigma2 = m*PSI+b, subject to constraints:
                    ##     {sigma2=Nt when PSI=0
                    ##      sigma2=10 when PSI=1}
                    ##   Then solve for m and b.
                    m = 10-Nt;
                    b=Nt;
                    sigma2 = m*PSI+b;
                    P2=-1;
                    while(P2<1 | P2>=maxP-1){ #truncated normal distribution
                        ##Since AI is unsorted, and we want to draw P2 s.t. AI(P2) ~ AI(P1),
                        ## we can draw a rank near the rank of P1, and convert back to index
                        P2r = round(rnorm(1,which(order(AI[1:(maxP-1)])==P1), sigma2)); #AI rank of P2
                        if(P2r>=1 & P2r<maxP-1){
                            P2 = order(AI[1:(maxP-1)])[P2r]; #transform to individual's index
                        }
                    }
                    
                    ## set AI of this individual as average of parents
                    AI[i] = (AI[P1]+AI[P2])/2;
                    if(is.na(AI[P1]) | is.na(AI[P2])){
                        cat("\n\nError\nit=",it,"; Nt=",Nt,"; maxP=",maxP,"; P1=",P1,": ",AI[P1],"; P2=",P2,": ",AI[P2],"\n");
                        return();
                    }
                }
            }
            N = Nt;
            M = min(M*(1+MG), 1);
            MAI = MAIt;
            AIout[sim, gen+1] = mean(AI,na.rm=T);
            ##cat("\n",sum(is.na(AI)))
        }
        setTxtProgressBar(pb, sim)
    }
    return(AIout);
}

## runBootstrap() runs bootstraps on the inferred AI data from Spear et al (2020) data.
runBootstrap = function(nBOOTS){# run bootstraps
    dat = read.table("../mexicans_sims_data.txt",header=T)
    ylim=c(0.39,0.55)*100;
    YL = ylim;
    plot(1,1,type='n',xlim=range(dat$Birth_year), ylim=ylim,xaxt='n',yaxt='n',xlab="",ylab="")
    axis(side=1,padj=-2, cex.axis=0.7)
    axis(side=2,padj=1.3,cex.axis=0.7)
    mtext(side=1,"Birth Year",line=1)
    mtext(side=2,"AI Ancestry (%)",line=1.15)
    
    ##include a status bar
    pb <- txtProgressBar(min = 0, max = nBOOTS, style = 3)
    
    ## bootstrap and add LOESS
    bd = array(NA, dim=c(2,nBOOTS,nrow(dat)));
    
    x01s = matrix(NA,nr=nBOOTS,nc=nrow(dat));
    p01s = matrix(NA,nr=nBOOTS,nc=nrow(dat));
    for(i in 1:nBOOTS){
        x=sample(1:nrow(dat),nrow(dat),replace=T)
        x01 = dat$Birth_year[x];
        x01s[i,] = x01;
        y01 = 100*dat$AI[x];
        SPAN=0.75
        l01 = loess(y01[order(x01)] ~ x01[order(x01)], span=SPAN)
        p01 = predict(l01, data.frame(x=x01[order(x01)],y=y01[order(x01)]));
        p01s[i,] = p01;
        lines(x01[order(x01)], p01, col="lightgrey");
        setTxtProgressBar(pb, i)
        yl = range(p01);
        if(YL[1] > yl[1]){
            YL[1] = yl[1];
        }
        if(YL[2] < yl[2]){
            YL[2] = yl[2];
        }
    }
    bd[1,,] = x01s;
    bd[2,,] = p01s;
    return(bd);
}

## this function takes parameters, runs the according simulation, and creates a stand-alone plot.
## plotAIsim() takes an additional parameter BGboot which will add background bootstraps of AI
##  inferred from the Spear et al (2020) paper. These bootstraps should be run separately and
##  included in the argument BGbootData.
plotAIsim <- function(NSIMS=100, NGENS=2, N0=1000, G=0, AI0=0.42, M0=0, MG=0, MAI0=0, MAIG=0,
                  FAI=0, PSI=0, BGboot=0, BGbootData=NA){

    if(BGboot == 1){
        if(is.null(dim(BGbootData))){
            cat("must run bootstraps with runBootstrap()\n");
            return();
        }
    
    
        ##rename bootstrap data and transform birth year to [0,2] generations
        bd = BGbootData;
        bd[1,,] = bd[1,,]-min(bd[1,,]);
        bd[1,,] = 2*bd[1,,]/max(bd[1,,]);
    }
    
    AIout = simAI(NSIMS=NSIMS, NGENS=NGENS, N0=N0, G=G, AI0=AI0, M0=M0, MG=MG, MAI0=MAI0,
                  MAIG=MAIG, FAI=FAI, PSI=PSI);
    
    FILEPREF=paste("AItraj_N0=",N0,"_G=",G,"_AI0=",AI0,"_M0=",M0,"_MG=",MG,"_MAI0=",MAI0,"_MAIG=",
                   MAIG,"_FAI=",FAI,"_PSI=",PSI,sep="");
    
    jpeg(paste(FILEPREF,".jpg",sep=""),height=2.5,width=2.5,pointsize=10,units="in",res=300)
    par(mar=c(2,2.2,0.3,0.3))
    plot(1,1,type='n',xlim=c(0,NGENS), ylim=c(0.37,0.6)*100,xaxt='n',yaxt='n',xlab="",ylab="")
    axis(side=1,padj=-2, cex.axis=0.7)
    axis(side=2,padj=1.3,cex.axis=0.7)
    mtext(side=1,"Generation",line=1)
    mtext(side=2,"AI anc (%)",line=1.15)
    
    if(BGboot == 1){
        nBOOTS = dim(BGbootData)[2];
        for(i in 1:nBOOTS){
            lines(bd[1,i,order(bd[1,i,])], bd[2,i,], col="lightgray");
        }
        legend("topleft",c("Bootstrap from data","Simulated"), col=c("lightgray",rgb(0,0,1,1)),
               lty=1,cex=0.8,box.col=rgb(1,1,1,0.8), bg=rgb(1,1,1,0.8));
    }
    
    apply(AIout*100,1,function(x){lines(0:NGENS, x, col=rgb(0,0,1,1))});

    ##add simulated means
    for(i in 0:NGENS){
        points(i, mean(AIout[,i+1])*100, pch=20,col='red')
    }
    box();
    dev.off();
    apply(AIout,2,mean)
    apply(AIout,2,sd)
}

                        
if(0){ #used for plotting random 1-off things
    ## SNM, nothing expected to happen
    plotAIsim(BGboot=1, BGbootData=bd);

    if(0){##plot population growth figures
        jpeg("ExpPopSize.jpg",height=2.5,width=3.75,pointsize=10,units="in",res=300)
        par(mar=c(2,2.2,0.3,0.3))
        plot(1,1,type='n',xlim=c(0,NGENS), ylim=c(N0,N0*(2)^2),xaxt='n',yaxt='n',xlab="",ylab="")
        axis(side=1,padj=-2, cex.axis=0.7)
        axis(side=2,padj=1.3,cex.axis=0.7)
        mtext(side=1,"Generation",line=1)
        mtext(side=2,"Population Size",line=1.15)
        cls = c("black","red","blue","grey");
        GR=c(0, 0.1, 0.5, 1);
        for(gi in 1:length(GR)){
            gr = GR[gi];
            n=N0;
            for(g in 1:2){
                n[g+1] = n[g]*(1+gr);
            }
            lines(0:2,n,col=cls[gi]);
            points(0:2,n,col=cls[gi],pch=20);
        }
        legend("topleft",rev(c("G=0","G=0.1","G=0.5","G=1")),col=rev(cls),lty=1);
        dev.off()

        Gset = c(0,0.1,0.5,1);
        for(G in Gset){
            plotAIsim(G=G, BGboot=1, BGbootData=bd);
        }
        
    }
    if(0){##plot assortative mating
        AI = sort(rbeta(N0, 4*AI0/(1-AI0), 4)) #initial population
        npar = 1000;
        PSI = c(0, 0.5, 0.75, 1);

        for(pi in 1:length(PSI)){
            p = PSI[pi];
            plotAIsim(PSI=p, BGboot=1, BGbootData=bd);

            AM = array(NA,dim=c(npar,2))
            for(i in 1:npar){
                P1 = sample(N0,1); #randomly choose parent 1
                m = 10-Nt;
                b=Nt;
                sigma2 = m*p+b;
                P2=-1;
                while(P2<1 | P2>N){ #truncated normal distribution
                    P2 = round(rnorm(1,P1, sigma2));
                }
                AM[i,] = c(AI[P1], AI[P2]);
            }

            jpeg(paste("AssortMat=",p,".jpg",sep=""),height=2.5,width=3.75,pointsize=10,units="in",res=300)
            par(mar=c(2,2.2,0.3,0.3))
            plot(1,1,type='n',xlim=range(AI), ylim=range(AI),xaxt='n',yaxt='n',xlab="",ylab="")
            axis(side=1,padj=-2, cex.axis=0.7)
            axis(side=2,padj=1.3,cex.axis=0.7)
            mtext(side=1,"AI(P1)",line=1)
            mtext(side=2,"AI(P2)",line=1.15)
            abline(0,1)
            points(AM[,1], AM[,2],col=rgb(0,0,0,0.2),pch=20);
            text(min(AM),max(AM),paste("cor=",format(cor(AM[,1],AM[,2]),digits=4),sep=""),pos=4)
            dev.off();
        }
    }

    if(0){ ## plot fecundity variance
        N=N0;
        AI = sort(rbeta(N, 4*AI0/(1-AI0), 4)) #initial population

        FAIr = c(0,0.1,0.2,0.4,0.8,1)
        jpeg(paste("fecundityDist.jpg",sep=""),height=2.5,width=3.75,pointsize=10,units="in",res=300)
        par(mar=c(2,2.2,0.3,0.3))
        yl=c(0,0.032);
        plot(1,1,type='n',xlim=c(1,N), ylim=yl,xaxt='n',yaxt='n',xlab="",ylab="")
        axis(side=1,padj=-2, cex.axis=0.7)
        axis(side=2,padj=1.3,cex.axis=0.7)
        mtext(side=1,"indiv (low-high AI)",line=1)
        mtext(side=2,"Prob Reprod",line=1.15)
        for(i in 1:length(FAIr)){
            FAI=FAIr[i];
            P1 = round((N-1)*rbeta(100000,1,1/(1+FAI)))+1; #draw P1 from (potentially) biased sample of domestic pop
            lines(density(P1),col=pcol[i])
        }
        legend("topleft",c("FAI=0","FAI=0.1","FAI=0.2","FAI=0.4","FAI=0.8","FAI=1.0"),col=pcol,lty=1)
        dev.off()
        for(f in FAIr){
            plotAIsim(FAI=f, PSI=0.75, BGboot=1, BGbootData=bd);
        }
    }

    if(0){ ## plot migration AI %
        
        jpeg(paste("migrantAI.jpg",sep=""),height=2.5,width=3.75,pointsize=10,units="in",res=300)
        par(mar=c(2,2.2,0.3,0.3))

        yl=c(0,3);
        plot(1,1,type='n',xlim=c(0,1), ylim=yl,xaxt='n',yaxt='n',xlab="",ylab="")
        axis(side=1,padj=-2, cex.axis=0.7)
        axis(side=2,padj=1.3,cex.axis=0.7)
        mtext(side=1,"AI",line=1)
        mtext(side=2,"Density",line=1.15)
        mai = c(0,0.1,0.2,0.4);
        for(i in 1:length(mai)){
            t = mai[i] + (1-mai[i])*AI0; #weighted average of current AI and 1.
            ait = rbeta(100000, 4*t/(1-t), 4)
            lines(density(ait),col=cls[i]);
        }
        legend("topleft",c("mAI=0","mAI=0.1","mAI=0.2","mAI=0.4"),col=cls,lty=1,cex=0.7)
        dev.off();

        Mset = c(0.1, 0.2, 0.4, 0.8);
        MAIset = c(0.1, 0.2, 0.4);
        FAIset = c(0, 0.1, 0.2, 0.4, 0.8, 1)
        AMset = c(0, 0.75)
        for(m in Mset){
            for(mai in MAIset){
                for(fai in FAIset){
                    for(am in AMset){
                        cat("\nM0=",m,"; MAI0=",mai,"; FAI=",fai,"; PSI=",am,"\n");
                        plotAIsim(M0=m, MAI0=mai, FAI=fai, PSI=am, BGboot=1, BGbootData=bd);
                    }
                }
            }
        }
    }

    if(0){ ## play with Melissa's AI data
        dat = read.table("../mexicans_sims_data.txt",header=T)
        jpeg(paste("AIdata.jpg",sep=""),height=2.5,width=3.75,pointsize=10,units="in",res=300)
        par(mar=c(2,2.2,0.3,0.3))
        plot(1,1,type='n',xlim=range(dat$Birth_year), ylim=range(100*dat$AI),xaxt='n',yaxt='n',xlab="",ylab="")
        plot(1,1,type='n',xlim=range(dat$Birth_year), ylim=range(100*c(0.38,0.55)),xaxt='n',yaxt='n',xlab="",ylab="")
        axis(side=1,padj=-2, cex.axis=0.7)
        axis(side=2,padj=1.3,cex.axis=0.7)
        mtext(side=1,"Birth Year",line=1)
        mtext(side=2,"AI Ancestry (%)",line=1.15)
        points(dat$Birth_year[which(dat$US_BORN==0)], 100*dat$AI[which(dat$US_BORN==0)], pch=16, col=rgb(1,0,1,0.1))
        points(dat$Birth_year[which(dat$US_BORN==1)], 100*dat$AI[which(dat$US_BORN==1)], pch=16, col=rgb(0,1,1,0.1))

        ##add regression lines for US born and not US born participants
        m01 = lm(100*dat$AI ~ dat$Birth_year)
        m0 = lm(100*dat$AI[which(dat$US_BORN==0)] ~ dat$Birth_year[which(dat$US_BORN==0)])
        m1 = lm(100*dat$AI[which(dat$US_BORN==1)] ~ dat$Birth_year[which(dat$US_BORN==1)])
        abline(m0, col=rgb(1,0,1,1),lwd=2)
        abline(m1, col=rgb(0,1,1,1),lwd=2)
        abline(m01, col=rgb(0,0,0,1),lwd=2)

        ## add loess curves for two groups
        x0 = dat$Birth_year[which(dat$US_BORN==0)];
        y0 = 100*dat$AI[which(dat$US_BORN==0)];
        x1 = dat$Birth_year[which(dat$US_BORN==1)];
        y1 = 100*dat$AI[which(dat$US_BORN==1)];
        x01 = dat$Birth_year;
        y01 = 100*dat$AI;
        SPAN=0.75
        l0 = loess(y0[order(x0)] ~ x0[order(x0)], span=SPAN)
        l1 = loess(y1[order(x1)] ~ x1[order(x1)], span=SPAN)
        l01 = loess(y01[order(x01)] ~ x01[order(x01)], span=0.5)
        p0 = predict(l0, data.frame(x=x0[order(x0)],y=y0[order(x0)]));
        p1 = predict(l1, data.frame(x=x1[order(x1)],y=y1[order(x1)]));
        p01 = predict(l01, data.frame(x=x01[order(x01)],y=y01[order(x01)]));
        lines(x0[order(x0)], p0, col=rgb(1,0,1,1),lwd=2)
        lines(x1[order(x1)], p1, col=rgb(0,1,1,1),lwd=2);
        lines(x01[order(x01)], p01, col=rgb(0,0,0,1),lwd=2);
        
        legend("bottomright",c("US born", "Not US born"), col=c(rgb(0,1,1,1), rgb(1,0,1,1)), pch=20, ncol=2, box.col=rgb(1,1,1,0.8), bg=rgb(1,1,1,0.8), cex=0.8)
        box()
        dev.off()
    }

    if(0){ ## play with Melissa's AI data - bootstrap loess curves
        nBOOTS=1000;
        bd = runBootstrap(nBOOTS);
        if(1){#plot it
            jpeg(paste("AIdata_bootstrapLOESS.jpg",sep=""),height=2.5,width=3.75,pointsize=10,units="in",res=300)
            par(mar=c(2,2.2,0.3,0.3))
            ylim=c(37,55);
            plot(1,1,type='n',xlim=range(dat$Birth_year), ylim=ylim,xaxt='n',yaxt='n',xlab="",ylab="")
            axis(side=1,padj=-2, cex.axis=0.7)
            axis(side=2,padj=1.3,cex.axis=0.7)
            mtext(side=1,"Birth Year",line=1)
            mtext(side=2,"AI Ancestry (%)",line=1.15)
            for(i in 1:nLOESS){
                lines(x01s[i,order(x01s[i,])], p01s[i,], col="lightgrey");
            }
            ##add regression lines
            m01 = lm(100*dat$AI ~ dat$Birth_year)
            abline(m01, col=rgb(1,1,1,1),lwd=2)
            abline(m01, col=rgb(0,0,1,1),lwd=2)
            if(YL[1] < ylim[1] | YL[2]>ylim[2]){
                cat("new ylim=c(",YL[1],",",YL[2],")\n");
            }
            box()
            dev.off()
        }
    }

        
}
