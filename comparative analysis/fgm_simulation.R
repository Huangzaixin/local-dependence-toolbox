##################################################################
###       Generate random data using the FGM copula model      ###
##################################################################
require(copula)

rndcopula <- rCopula(300000, fgmCopula(0.9, dim=2))
u <- rndcopula[,1]
v <- rndcopula[,2]

# save data
write.csv(u, file = "fgm_u.csv", row.names = FALSE)
write.csv(v, file = "fgm_v.csv", row.names = FALSE)

