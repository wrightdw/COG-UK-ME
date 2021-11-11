install.packages("twitteR")
install.packages(c("devtools", "rjson", "bit64", "httr"))
install.packages("base64enc")
devtools::install_github("jrowen/twitteR", ref = "oauth_httr_1_0")
library(base64enc)
library(devtools)
library(rjson)
library(bit64)
library(httr)
library(twitteR)

consumer_key <- "dVjqTfZ75Y8GiqWXlfRY4Gdi3"
access_token<- "1408421851408355332-ijmf06Vetwpr7H0agk4LQ4YGlLOpbv"
consumer_secret <- "RWjXDoOdFpT2T2MpGHgShgW621GrPNiRx711PUdbmG0UOjDbTo"
access_secret<- "G5TIZNtOD1aRRaVuXy51FMIiag8oOysIEhiHQxtMbjjIy"
setup_twitter_oauth(consumer_key, consumer_secret,access_token,access_secret)

tweet("iphone")
tw<-updateStatus('We will publish new plots very soon')
## Ronapreve tweet
# generate png file

png(filename="Ronapreve.png",width=1200, height=650)
ronapreve_upset
dev.off()
tweet("Monitoring the emergence of new mutations in the UK that might affect the mAb cocktail. For more information please have a look at http://sars2.cvr.gla.ac.uk/cog-uk/", mediaPath = "Ronapreve.png")


