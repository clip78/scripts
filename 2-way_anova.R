with(data, tapply(score, list(FactorA,FactorB), mean))

hill$test = factor(hill$test, levels=c("pre","post"))

with(hill, tapply(SSS, list(diet,test), mean))

> with(hill, boxplot(SSS ~ diet + test))    # output not shown
> with(hill, boxplot(SSS ~ test + diet))    # compare this one with the last one
> title(main="Ben Hill's SSS Data")
> title(ylab="SSS Scores")

> aov.out = aov(SSS ~ diet * test + Error(subject/test), data=hill)
> summary(aov.out)

Other Designs

For reference, here are model formulae for a couple other common designs... 

Two factor design with repeated measures on both factors:
DV ~ IV1 * IV2 + Error(subject/(IV1*IV2))

Three factor mixed design with repeated measures on IV2 and IV3:
DV ~ IV1 * IV2 * IV3 + Error(subject/(IV2*IV3))