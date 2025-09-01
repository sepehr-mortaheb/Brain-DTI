library(lme4)
library(emmeans)
library(magrittr)
library(Matrix)
library(car)
library(lmerTest)

## Whole-Brain FC-SC Corrlation Analysis 
fcsc_dir = '/Users/sepehrmortaheb/Desktop/Results/Results 2/SCFC_Analysis/'
sim_dir = paste(fcsc_dir, 'scfc_corr_parrcorr_HO.csv', sep="")
df = read.csv(sim_dir)

df$subject = as.factor(df$subject)
df$group = as.factor(df$group)
df$time_point = as.factor(df$time_point)
df$flight = as.factor(df$flight)
filt_df = df[(df$flight=='f1') & (df$time_point!='tp3'),]

mix = lmer( corr ~ group*time_point + (1|subject), data=filt_df)
Anova(mix)
mix.emm = emmeans(mix, "time_point", by="group")
mix.pairs = pairs(mix.emm)
mix.pairs

fxm = lm( corr ~ group*time_point, data=filt_df)
Anova(fxm)
fxm.emm = emmeans(fxm, "time_point", by="group")
fxm.pairs = pairs(fxm.emm)
fxm.pairs
############################################################################
## Whole-Brain FC-SC Cosine Similarity Analysis 
fcsc_dir = '/Users/sepehrmortaheb/Desktop/Results/Results 2/SCFC_Analysis/'
sim_dir = paste(fcsc_dir, 'scfc_corr_parrcorr.csv', sep="")
df = read.csv(sim_dir)

df$subject = as.factor(df$subject)
df$group = as.factor(df$group)
df$time_point = as.factor(df$time_point)
df$flight = as.factor(df$flight)
filt_df = df[(df$flight=='f1') & (df$time_point!='tp3'),]

mix = lmer( cosine_sim ~ group*time_point + (1|subject), data=filt_df)
Anova(mix)
mix.emm = emmeans(mix, "time_point", by="group")
mix.pairs = pairs(mix.emm)
mix.pairs

fxm = lm( cosine_sim ~ group*time_point, data=filt_df)
Anova(fxm)
fxm.emm = emmeans(fxm, "time_point", by="group")
fxm.pairs = pairs(fxm.emm)
fxm.pairs
##########################################################################
## Whole-Brain Liberality Analysis 
gsp_dir = '/Users/sepehrmortaheb/Desktop/Results/Results 2/GSP_Analysis/'
alglib_dir = paste(gsp_dir, 'df_Lib_Alg_whole_HO.csv', sep="")
df = read.csv(alglib_dir)

df$subject = as.factor(df$subject)
df$group = as.factor(df$group)
df$time_point = as.factor(df$time_point)
df$flight = as.factor(df$flight)
filt_df = df[(df$flight==1) & (df$time_point!=3),]

mix = lmer( Liberality ~ group*time_point + (1|subject), data=filt_df)
Anova(mix)
mix.emm = emmeans(mix, "time_point", by="group")
mix.pairs = pairs(mix.emm)
mix.pairs

fxm = lm( Liberality ~ group*time_point , data=filt_df)
Anova(fxm)
fxm.emm = emmeans(fxm, "time_point", by="group")
fxm.pairs = pairs(fxm.emm)
fxm.pairs
##########################################################################
## Whole-Brain Alignment Analysis 
gsp_dir = '/Users/sepehrmortaheb/Desktop/Results/Results 2/GSP_Analysis/'
alglib_dir = paste(gsp_dir, 'df_Lib_Alg_whole_HO.csv', sep="")
df = read.csv(alglib_dir)

df$subject = as.factor(df$subject)
df$group = as.factor(df$group)
df$time_point = as.factor(df$time_point)
df$flight = as.factor(df$flight)
filt_df = df[(df$flight==1) & (df$time_point!=3),]

mix = lmer( Alignmnet ~ group*time_point + (1|subject), data=filt_df)
Anova(mix)
mix.emm = emmeans(mix, "time_point", by="group")
mix.pairs = pairs(mix.emm)
mix.pairs

fxm = lm( Alignmnet ~ group*time_point , data=filt_df)
Anova(fxm)
fxm.emm = emmeans(fxm, "time_point", by="group")
fxm.pairs = pairs(fxm.emm)
fxm.pairs
####################################################################
## Regional Coupling Analysis 
gsp_dir = '/Users/sepehrmortaheb/Desktop/Results/Results 2/GSP_Analysis/'
file_dir = paste(gsp_dir, 'df_SDI_HO.csv', sep="")
df = read.csv(file_dir)

df$subject = as.factor(df$subject)
df$group = as.factor(df$group)
df$time_point = as.factor(df$time_point)
df$flight = as.factor(df$flight)
filt_df = df[(df$flight==1) & (df$time_point!=3) & (df$region==52),]

mix = lmer( c_index ~ group*time_point + (1|subject), data=filt_df)
Anova(mix)
mix.emm = emmeans(mix, "time_point", by="group")
mix.pairs = pairs(mix.emm)
mix.pairs

fxm = lm( c_index ~ group*time_point , data=filt_df)
Anova(fxm)
fxm.emm = emmeans(fxm, "time_point", by="group")
fxm.pairs = pairs(fxm.emm)
fxm.pairs
#####################################################################
## Regional Deoupling Analysis 
gsp_dir = '/Users/sepehrmortaheb/Desktop/Results/Results 2/GSP_Analysis/'
file_dir = paste(gsp_dir, 'df_SDI_HO.csv', sep="")
df = read.csv(file_dir)

df$subject = as.factor(df$subject)
df$group = as.factor(df$group)
df$time_point = as.factor(df$time_point)
df$flight = as.factor(df$flight)
filt_df = df[(df$flight==1) & (df$time_point!=3) & (df$region==52),]

mix = lmer( d_index ~ group*time_point + (1|subject), data=filt_df)
Anova(mix)
mix.emm = emmeans(mix, "time_point", by="group")
mix.pairs = pairs(mix.emm)
mix.pairs

fxm = lm( d_index ~ group*time_point , data=filt_df)
Anova(fxm)
fxm.emm = emmeans(fxm, "time_point", by="group")
fxm.pairs = pairs(fxm.emm)
fxm.pairs
#####################################################################
## Regional SDI Analysis 
gsp_dir = '/Users/sepehrmortaheb/Desktop/Results/Results 2/GSP_Analysis/'
file_dir = paste(gsp_dir, 'df_SDI_HO.csv', sep="")
df = read.csv(file_dir)

df$subject = as.factor(df$subject)
df$group = as.factor(df$group)
df$time_point = as.factor(df$time_point)
df$flight = as.factor(df$flight)
filt_df = df[(df$flight==1) & (df$time_point!=3) & (df$region==52),]

mix = lmer( sdi ~ group*time_point + (1|subject), data=filt_df)
Anova(mix)
mix.emm = emmeans(mix, "time_point", by="group")
mix.pairs = pairs(mix.emm)
mix.pairs

fxm = lm( sdi ~ group*time_point , data=filt_df)
Anova(fxm)
fxm.emm = emmeans(fxm, "time_point", by="group")
fxm.pairs = pairs(fxm.emm)
fxm.pairs

#################################################################
## All Regions Coupling Analysis 
gsp_dir = '/Users/sepehrmortaheb/Desktop/Results/Results 2/GSP_Analysis/'
file_dir = paste(gsp_dir, 'df_SDI_HO.csv', sep="")
df = read.csv(file_dir)

df$subject = as.factor(df$subject)
df$group = as.factor(df$group)
df$time_point = as.factor(df$time_point)
df$flight = as.factor(df$flight)

for (i in 0:62){
  filt_df = df[(df$flight==1) & (df$time_point!=3) & (df$region==i),]
  mix = lmer( c_index ~ group*time_point + (1|subject), data=filt_df)
  res = Anova(mix)
  if (res[,"Pr(>Chisq)"][3]<0.05){
    print(i)
    print(res)
  }
}

#################################################################
## All Regions Decoupling Analysis 
gsp_dir = '/Users/sepehrmortaheb/Desktop/Results/Results 2/GSP_Analysis/'
file_dir = paste(gsp_dir, 'df_SDI_HO.csv', sep="")
df = read.csv(file_dir)

df$subject = as.factor(df$subject)
df$group = as.factor(df$group)
df$time_point = as.factor(df$time_point)
df$flight = as.factor(df$flight)

for (i in 0:62){
  filt_df = df[(df$flight==1) & (df$time_point!=3) & (df$region==i),]
  mix = lmer( d_index ~ group*time_point + (1|subject), data=filt_df)
  res = Anova(mix)
  if (res[,"Pr(>Chisq)"][3]<0.05){
    print(i)
    print(res)
  }
}

#################################################################
## All Regions sdi Analysis 
gsp_dir = '/Users/sepehrmortaheb/Desktop/Results/Results 2/GSP_Analysis/'
file_dir = paste(gsp_dir, 'df_SDI_HO.csv', sep="")
df = read.csv(file_dir)

df$subject = as.factor(df$subject)
df$group = as.factor(df$group)
df$time_point = as.factor(df$time_point)
df$flight = as.factor(df$flight)

for (i in 0:62){
  filt_df = df[(df$flight==1) & (df$time_point!=3) & (df$region==i),]
  mix = lmer( sdi ~ group*time_point + (1|subject), data=filt_df)
  res = Anova(mix)
  if (res[,"Pr(>Chisq)"][3]<0.05){
    print(i)
    print(res)
  }
}