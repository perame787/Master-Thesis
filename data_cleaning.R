library(quantmod)
library(data.table)
tickers<-c("IVV", "EZU", "JPXN", "EPP", "CEC.PA", "EEM", "VCIT", "EMB", "TRET.L", "GSG")
# SG CTA Index: NEIXCTA
getSymbols(Symbols = tickers, from="2011-06-01", index.return = TRUE)

#Convert xts objects to data table and add Date as column
for(i in seq_along(tickers)){
    assign(tickers[i], as.data.table(get(tickers[i])))
}

#Load time series that were downloaded sep. from Bloomberg
DE10Y<-fread("C:/Users/Pedro/Documents/Studium/WU Dokummente/SoSe 23/Thesis/Scripts/Data/DE_FI.csv")
NEIXCTAT<-fread("C:/Users/Pedro/Documents/Studium/WU Dokummente/SoSe 23/Thesis/Scripts/Data/NEIXCTAT.csv")
USFI<-fread("C:/Users/Pedro/Documents/Studium/WU Dokummente/SoSe 23/Thesis/Scripts/Data/US_FI.csv")

#Compute bond prices for FI instruments
USFI10Y<-USFI[,c("Date", "USGG10YR")]
USFI10Y$USGG10YR<-USFI10Y$USGG10YR/100
USFI10Y$USGG10YR<-100/(1+USFI10Y$USGG10YR)^10

DE10Y$GTDEM10YR<-DE10Y$GTDEM10YR/100
DE10Y$GTDEM10YR<-100/(1+DE10Y$GTDEM10YR)^10

bloomberg<-list(DE10Y, USFI10Y, NEIXCTAT)

#Build data set with just adj. close prices of assets
prices_ts<-function(tickers, as.xts=TRUE){
  if(as.xts){
    prices<-na.omit(merge(get(tickers[1])[, 4], 
                          get(tickers[2])[, 4]))
    for(i in 3:(length(tickers))){
      prices<-na.omit(merge(prices, get(tickers[i])[, 4]))
    } 
  }
  else{
    prices<-na.omit(merge(get(tickers[1])[, c(1, 7)], 
                          get(tickers[2])[, c(1, 7)], 
                          by=c("index"), all=T))
    for(i in 3:(length(tickers))){
      prices<-na.omit(merge(prices, get(tickers[i])[, c(1,7)], by=c("index"), all=T))
    } 
  }
  colnames(prices)[1]<-"Date"
  prices
}
#undebug(prices_ts)
#Prices of the scenario set for TS from Yahoo
prices_ss<-prices_ts(tickers, as.xts = F)

#Add TS from Bloomberg
for(i in bloomberg){
  i$Date<-as.Date(i$Date)
  prices_ss<-na.omit(merge(prices_ss, i, by = "Date", all = T))
}
tickers<-c(tickers, "DE10Y", "USFI10Y", "NEIXCTAT")

#Build empirical scenario set using log returns*(-1) of the prices_ss df
loss_dist<-function(price_df, names, ret_type=c("log", "simple"), no_zeroes=TRUE){
  #Number of rows - 1 because of NAs when computing returns
  df<-as.data.table(matrix(NA, nrow = (nrow(price_df)-1), ncol = ncol(price_df)))
  colnames(df)<-c("Date", names)
  df[,1]<-price_df[,1][-1]
  for(i in 2:ncol(df)){
    if(ret_type=="log"){
      df[,i]<-na.omit(diff(log(unlist(price_df[, i, with = F]))))*(-1)
    }
    else if(ret_type=="simple"){
      df[,i]<-(exp(na.omit(diff(log(unlist(price_df[, i, with = F])))))-1)*(-1) 
    }
  }
  df<-df[, -1]
  if(no_zeroes==TRUE){
    df<-df[apply(df!=0, 1, all),]
  }
  df
}

#Loss distributions-Alternative function due to necessity of xts object
#loss_dist_xts<-function(price_df, names, ret_type=c("log", "simple")){
for(i in 1:ncol(price_df)){
    if(ret_type=="log"){
      price_df[,i]<-c(0, na.omit(diff(log(unlist(as.vector(price_df[, i, with = F])))))*(-1))
    }
    else{
      price_df[,i]<-c(0, (exp(na.omit(diff(log(unlist(as.vector(price_df[, i, with = F]))))))-1)*(-1))
    }
  }
price_df

#undebug(loss_dist)
scenario_set_log<-loss_dist(price_df=prices_ss, names=tickers, ret_type = "log", no_zeroes = T)
scenario_set_simple<-loss_dist(price_df=prices_ss, names=tickers, ret_type = "simple", no_zeroes = T)
#Remove 
names(prices_ss)<-c("Date", tickers)
write.csv(prices_ss, "C:/Users/Pedro/Documents/Studium/WU Dokummente/SoSe 23/Thesis/Scripts/Data/prices_ss.csv", row.names = F)
write.csv(scenario_set_log, "C:/Users/Pedro/Documents/Studium/WU Dokummente/SoSe 23/Thesis/Scripts/Data/scenario_set_log.csv", row.names = F)
write.csv(scenario_set_simple, "C:/Users/Pedro/Documents/Studium/WU Dokummente/SoSe 23/Thesis/Scripts/Data/scenario_set_simple.csv", row.names = F)
