# test only works with package drcSeedGerm

#One single germination curve ##########################################
# library(drcSeedGerm)
# mod <- drm(count ~ start + end, data=chickweed, type="event", fct=LL.3())
# summary(mod)
# plot(mod)
# 
# #Multiple germination curves ##########################################
# data(verbascum)
# # mod <- drm(nSeeds ~ timeBef + timeAf, data=verbascum, fct=LL.3(),
# #              curveid=Species, type="event" )
# 
# #The error is due to the fact that the self starting routine returns a 'd' 
# #value for the third species that is higher than 1.
# #The error can be overriden by manually supplying correct starting values
# verbascum_cum <- subset(verbascum, is.finite(timeAf)==T)
# mod <- drm(I(nCum/25) ~ timeAf, data=verbascum, fct=LL.3(),
#            curveid=Species)
# mod2 <- drm(nSeeds ~ timeBef + timeAf, data=verbascum, fct=LL.3(),
#             curveid=Species, type="event", start=coef(mod))
# summary(mod2)
# plot(mod2)
# 
# #Hydrotime-to-event model
# data(rape)
# HTEmod <- drm(nSeeds ~ timeBef + timeAf + Psi, data = rape, fct=HTE1(), type = "event")
# summary(HTEmod)
# 
# #Thermal-time-to-event model
# data(barley)
# TTEmod <- drm(nSeeds ~ timeBef + timeAf + Temp, 
#                    fct=TTERF(), data = barley, type="event")
# summary(TTEmod)

