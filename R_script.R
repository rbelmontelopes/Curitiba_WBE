#R Script for the analysis performed by Belmonte-Lopes et al. "Wastewater-based epidemiology of SARS-COV-2 anticipated outbreak of COVID-19 cases in Curitiba (Brazil) before the first records of Omicron infections" submited as a letter to The Lancet
# read the wastewater data
data = read_csv("wastewater_data.csv")

data["N1 copies L"] = (data$`N1 copies per microL`*100)/data$`filtered volume(L)`
data["sd N1 copies L"] = (data$`sd N1 copies microL`*100)/data$`filtered volume(L)`

#calculate upper and lower 95% Confidence intervals
data["upper95CI_copiesL"] = data$`N1 copies L` + (1.96)*(data$`sd N1 copies L`)/sqrt(3)
data["lower95CI_copiesL"] = data$`N1 copies L` - (1.96)*(data$`sd N1 copies L`)/sqrt(3)

#calculate daily viral load and 95% CI
data["daily_viral_load"] = data$`N1 copies L` * (((data$flow*60)*60)*24)
data["upper95CI_daily_viral_load"] = data$upper95CI_copiesL* (((data$flow*60)*60)*24)
data["lower95CI_daily_viral_load"] = data$lower95CI_copiesL* (((data$flow*60)*60)*24)

# transform data into a time series
data = data %>% as_tsibble(key = WWTP)
#generate rows for missing days
data_int = fill_gaps(data)

#interpolate data for ETE-AT
for(i in c(3:13)){data_int[data_int[,2] == "ETE-AT",colnames(data_int[i])]=na.interp(data_int[data_int[,2]=="ETE-AT",][colnames(data_int[i])],lambda="auto")}
#interpolate data for ETE-BL
for(i in c(3:13)){data_int[data_int[,2] == "ETE-BL",colnames(data_int[i])]=na.interp(data_int[data_int[,2]=="ETE-BL",][colnames(data_int[i])],lambda="auto")}
#interpolate data for ETE-CX
for(i in c(3:13)){data_int[data_int[,2] == "ETE-CX",colnames(data_int[i])]=na.interp(data_int[data_int[,2]=="ETE-CX",][colnames(data_int[i])],lambda="auto")}
#interpolate data for ETE-PA
for(i in c(3:13)){data_int[data_int[,2] == "ETE-PA",colnames(data_int[i])]=na.interp(data_int[data_int[,2]=="ETE-PA",][colnames(data_int[i])],lambda="auto")}
#interpolate data for ETE-SQ
for(i in c(3:13)){data_int[data_int[,2] == "ETE-SQ",colnames(data_int[i])]=na.interp(data_int[data_int[,2]=="ETE-SQ",][colnames(data_int[i])],lambda="auto")}

# generate filters to get the data always from Tuesdays
Day = lubridate::wday(data_int$Date,label=TRUE)
Week = Day %in% c("Ter")
#get the data from thursdays
data_int = data_int[Week,]
#set the tsibble to 7D again
data_int = as_tsibble(data_int, index=Date)

#summarise the data by week
weekly = data_int %>%  index_by(Date) %>% 
  summarise(daily_viral_load = sum(daily_viral_load),
            upper95CI_daily_viral_load = sum(upper95CI_daily_viral_load),
            lower95CI_daily_viral_load = sum(lower95CI_daily_viral_load))


#calculate the log10 of the viral load and of its 95%CI
weekly["log10_daily_viral_load"] = log10(weekly$daily_viral_load)
weekly["log10_upper95CI_daily_viral_load"] = log10(weekly$upper95CI_daily_viral_load)
weekly["log10_lower95CI_daily_viral_load"] = log10(weekly$lower95CI_daily_viral_load)

#read the active cases data            
active_cases = read_csv("active_cases.csv")

#bind the weekly and active cases data
weekly = bind_cols(weekly,active_cases[,2:22])

#get the numerical data for correlation analysis
cor_data = weekly[,c(2,5,8:28)]

#rename columns for ploting
colnames(cor_data) = c("N1 copies/day", "log10(N1 copies/day)", "active cases", "active cases d+1", "active cases d+2", "active cases d+3", "active cases d+4", "active cases d+5", "active cases d+6", "active cases d+7", "active cases d+8", "active cases d+9", "active cases d+10", "active cases d+11", "active cases d+12", "active cases d+13", "active cases d+14", "active cases d+15", "active cases d+16", "active cases d+17", "active cases d+18", "active cases d+19", "active cases d+20", "active cases d+21")

#generate a correlation matrix
cor_mat_weekly = cor_mat(cor_data)
#plot the results
cor_plot(cor_mat_weekly, label=TRUE,font.label = list(size=0.6, color="white"))


# calculate the axis scale factor for the plot of active cases vs log10 viral load
#get the ranges of the data
r_ca = range(weekly$active_cases)
r_l10 = range(weekly$log10_daily_viral_load)
r_l10_U = range(weekly$log10_upper95CI_daily_viral_load)
r_l10_L = range(weekly$log10_lower95CI_daily_viral_load)

# calculate the scale factor
sf = diff(r_ca)/diff(r_l10) - 1700

#sf_U = diff(r_ca)/diff(r_l10_U)
#sf_L = diff(r_ca)/diff(r_l10_L)

# save the transformation to be used in an object
trans <- ~ ((. - r_ca[1]) / sf) + r_l10[1]

# generate some objects to use in plotting the text
detection = c("Delta detection", "Omicron detection")
sec_axs_comp = "viral load (N1 copies/day)"

#Generate the plot using ggplot
ggplot(weekly)+geom_rect(data=NULL,aes(xmin=as.Date("2021-03-14"),xmax=as.Date("2021-04-04"),ymin=-Inf,ymax=Inf),fill="red")+geom_rect(data=NULL,aes(xmin=as.Date("2021-05-29"),xmax=as.Date("2021-06-08"),ymin=-Inf,ymax=Inf),fill="red")+geom_rect(data=NULL,aes(xmin=as.Date("2021-03-02"),xmax=as.Date("2021-03-13"),ymin=-Inf,ymax=Inf),fill="#FFCC39")+geom_rect(data=NULL,aes(xmin=as.Date("2021-04-05"),xmax=as.Date("2021-05-28"),ymin=-Inf,ymax=Inf),fill="#FFCC39")+geom_rect(data=NULL,aes(xmin=as.Date("2021-06-09"),xmax=as.Date("2021-07-07"),ymin=-Inf,ymax=Inf),fill="#FFCC39")+geom_rect(data=NULL,aes(xmin=as.Date("2021-07-08"),xmax=as.Date("2022-02-24"), ymin=-Inf,ymax=Inf),fill="#FFFF99")+geom_vline(xintercept = as.Date("2021-03-02"), color="black", linetype="dashed", size=0.5) +annotate(geom="text", x=as.Date("2021-03-20"), y=16000, label="50% Gamma",color="black", size=3)+geom_vline(xintercept = as.Date("2021-04-01"), color="black", linetype="dashed", size=0.5) +annotate(geom="text", x=as.Date("2021-04-20"), y=14000, label="100% Gamma",color="black", size=3)+geom_vline(xintercept = as.Date("2021-05-13"), color="black", linetype="twodash", size=0.5) +annotate(geom="text", x=as.Date("2021-06-07"), y=16000, label=bquote(1^a ~ .(detection[1])),color="black", size=3)+geom_vline(xintercept = as.Date("2021-08-01"), color="black", linetype="twodash", size=0.5) +annotate(geom="text", x=as.Date("2021-08-15"), y=16000, label="45% Delta",color="black", size=3)+ geom_vline(xintercept = as.Date("2021-12-14"), color="black", linetype="dotted", size=0.5)+annotate(geom="text", x=as.Date("2022-01-11"), y=18000, label=bquote(1^a ~ .(detection[2])),color="black", size=3)+geom_col(aes(x=`Date`, y=`active_cases`), color="black", fill="white")+geom_ribbon(aes(x=`Date`, ymin=((`log10_lower95CI_daily_viral_load` - r_l10[1])*sf+r_ca[1]), ymax=((log10_upper95CI_daily_viral_load -r_l10[1])*sf+r_ca[1]), color="#CCCCCC", alpha=0.2 ), color="grey")+geom_line(aes(x=`Date`, y=((log10_daily_viral_load-r_l10[1])*sf+r_ca[1])), size=0.8, color="black")+scale_y_continuous(name = "COVID-19 Active cases", sec.axis = sec_axis(trans=trans, name = bquote(log[10] ~ .(sec_axs_comp[1]))))+scale_x_date(name="Date", breaks = seq(weekly$Date[1],weekly$Date[52], by=14), labels= format(seq(weekly$Date[1],weekly$Date[52], by=14), "%m-%d-%Y"))+theme(axis.text.x = element_text(angle = 45, hjust = 1, color="black"), axis.line = element_line(colour = "black"), panel.background = element_rect(fill = "transparent", colour = NA), legend.position = "none")


# read the wastewater data
data = read_csv("wastewater_data.csv")

data["N1 copies L"] = (data$`N1 copies per microL`*100)/data$`filtered volume(L)`
data["sd N1 copies L"] = (data$`sd N1 copies microL`*100)/data$`filtered volume(L)`

#calculate upper and lower 95% Confidence intervals
data["upper95CI_copiesL"] = data$`N1 copies L` + (1.96)*(data$`sd N1 copies L`)/sqrt(3)
data["lower95CI_copiesL"] = data$`N1 copies L` - (1.96)*(data$`sd N1 copies L`)/sqrt(3)

#calculate daily viral load and 95% CI
data["daily_viral_load"] = data$`N1 copies L` * (((data$flow*60)*60)*24)
data["upper95CI_daily_viral_load"] = data$upper95CI_copiesL* (((data$flow*60)*60)*24)
data["lower95CI_daily_viral_load"] = data$lower95CI_copiesL* (((data$flow*60)*60)*24)

# transform data into a time series
data = data %>% as_tsibble(key = WWTP)
#generate rows for missing days
data_int = fill_gaps(data)

#interpolate data for ETE-AT
for(i in c(3:13)){data_int[data_int[,2] == "ETE-AT",colnames(data_int[i])]=na.interp(data_int[data_int[,2]=="ETE-AT",][colnames(data_int[i])],lambda="auto")}
#interpolate data for ETE-BL
for(i in c(3:13)){data_int[data_int[,2] == "ETE-BL",colnames(data_int[i])]=na.interp(data_int[data_int[,2]=="ETE-BL",][colnames(data_int[i])],lambda="auto")}
#interpolate data for ETE-CX
for(i in c(3:13)){data_int[data_int[,2] == "ETE-CX",colnames(data_int[i])]=na.interp(data_int[data_int[,2]=="ETE-CX",][colnames(data_int[i])],lambda="auto")}
#interpolate data for ETE-PA
for(i in c(3:13)){data_int[data_int[,2] == "ETE-PA",colnames(data_int[i])]=na.interp(data_int[data_int[,2]=="ETE-PA",][colnames(data_int[i])],lambda="auto")}
#interpolate data for ETE-SQ
for(i in c(3:13)){data_int[data_int[,2] == "ETE-SQ",colnames(data_int[i])]=na.interp(data_int[data_int[,2]=="ETE-SQ",][colnames(data_int[i])],lambda="auto")}

# generate filters to get the data always from Tuesdays
Day = lubridate::wday(data_int$Date,label=TRUE)
Week = Day %in% c("Ter")
#get the data from thursdays
data_int = data_int[Week,]
#set the tsibble to 7D again
data_int = as_tsibble(data_int, index=Date)

#summarise the data by week
weekly = data_int %>%  index_by(Date) %>% 
  summarise(daily_viral_load = sum(daily_viral_load),
            upper95CI_daily_viral_load = sum(upper95CI_daily_viral_load),
            lower95CI_daily_viral_load = sum(lower95CI_daily_viral_load))


#calculate the log10 of the viral load and of its 95%CI
weekly["log10_daily_viral_load"] = log10(weekly$daily_viral_load)
weekly["log10_upper95CI_daily_viral_load"] = log10(weekly$upper95CI_daily_viral_load)
weekly["log10_lower95CI_daily_viral_load"] = log10(weekly$lower95CI_daily_viral_load)

#read the active cases data            
active_cases = read_csv("active_cases.csv")

#bind the weekly and active cases data
weekly = bind_cols(weekly,active_cases[,2:22])

#get the numerical data for correlation analysis
cor_data = weekly[,c(2,5,8:28)]

#rename columns for ploting
colnames(cor_data) = c("N1 copies/day", "log10(N1 copies/day)", "active cases", "active cases d+1", "active cases d+2", "active cases d+3", "active cases d+4", "active cases d+5", "active cases d+6", "active cases d+7", "active cases d+8", "active cases d+9", "active cases d+10", "active cases d+11", "active cases d+12", "active cases d+13", "active cases d+14", "active cases d+15", "active cases d+16", "active cases d+17", "active cases d+18", "active cases d+19", "active cases d+20", "active cases d+21")

#generate a correlation matrix
cor_mat_weekly = cor_mat(cor_data)
#plot the results
cor_plot(cor_mat_weekly, label=TRUE,font.label = list(size=0.6, color="white"))


# calculate the axis scale factor for the plot of active cases vs log10 viral load
#get the ranges of the data
r_ca = range(weekly$active_cases)
r_l10 = range(weekly$log10_daily_viral_load)
r_l10_U = range(weekly$log10_upper95CI_daily_viral_load)
r_l10_L = range(weekly$log10_lower95CI_daily_viral_load)

# calculate the scale factor
sf = diff(r_ca)/diff(r_l10) - 1700

#sf_U = diff(r_ca)/diff(r_l10_U)
#sf_L = diff(r_ca)/diff(r_l10_L)

# save the transformation to be used in an object
trans <- ~ ((. - r_ca[1]) / sf) + r_l10[1]

# generate some objects to use in plotting the text
detection = c("Delta detection", "Omicron detection")
sec_axs_comp = "viral load (N1 copies/day)"

#Generate the plot using ggplot
ggplot(weekly)+geom_rect(data=NULL,aes(xmin=as.Date("2021-03-14"),xmax=as.Date("2021-04-04"),ymin=-Inf,ymax=Inf),fill="red")+geom_rect(data=NULL,aes(xmin=as.Date("2021-05-29"),xmax=as.Date("2021-06-08"),ymin=-Inf,ymax=Inf),fill="red")+geom_rect(data=NULL,aes(xmin=as.Date("2021-03-02"),xmax=as.Date("2021-03-13"),ymin=-Inf,ymax=Inf),fill="#FFCC39")+geom_rect(data=NULL,aes(xmin=as.Date("2021-04-05"),xmax=as.Date("2021-05-28"),ymin=-Inf,ymax=Inf),fill="#FFCC39")+geom_rect(data=NULL,aes(xmin=as.Date("2021-06-09"),xmax=as.Date("2021-07-07"),ymin=-Inf,ymax=Inf),fill="#FFCC39")+geom_rect(data=NULL,aes(xmin=as.Date("2021-07-08"),xmax=as.Date("2022-02-24"), ymin=-Inf,ymax=Inf),fill="#FFFF99")+geom_vline(xintercept = as.Date("2021-03-02"), color="black", linetype="dashed", size=0.5) +annotate(geom="text", x=as.Date("2021-03-20"), y=16000, label="50% Gamma",color="black", size=3)+geom_vline(xintercept = as.Date("2021-04-01"), color="black", linetype="dashed", size=0.5) +annotate(geom="text", x=as.Date("2021-04-20"), y=14000, label="100% Gamma",color="black", size=3)+geom_vline(xintercept = as.Date("2021-05-13"), color="black", linetype="twodash", size=0.5) +annotate(geom="text", x=as.Date("2021-06-07"), y=16000, label=bquote(1^a ~ .(detection[1])),color="black", size=3)+geom_vline(xintercept = as.Date("2021-08-01"), color="black", linetype="twodash", size=0.5) +annotate(geom="text", x=as.Date("2021-08-15"), y=16000, label="45% Delta",color="black", size=3)+ geom_vline(xintercept = as.Date("2021-12-14"), color="black", linetype="dotted", size=0.5)+annotate(geom="text", x=as.Date("2022-01-11"), y=18000, label=bquote(1^a ~ .(detection[2])),color="black", size=3)+geom_col(aes(x=`Date`, y=`active_cases`), color="black", fill="white")+geom_ribbon(aes(x=`Date`, ymin=((`log10_lower95CI_daily_viral_load` - r_l10[1])*sf+r_ca[1]), ymax=((log10_upper95CI_daily_viral_load -r_l10[1])*sf+r_ca[1]), color="#CCCCCC", alpha=0.2 ), color="grey")+geom_line(aes(x=`Date`, y=((log10_daily_viral_load-r_l10[1])*sf+r_ca[1])), size=0.8, color="black")+scale_y_continuous(name = "COVID-19 Active cases", sec.axis = sec_axis(trans=trans, name = bquote(log[10] ~ .(sec_axs_comp[1]))))+scale_x_date(name="Date", breaks = seq(weekly$Date[1],weekly$Date[52], by=14), labels= format(seq(weekly$Date[1],weekly$Date[52], by=14), "%m-%d-%Y"))+theme(axis.text.x = element_text(angle = 45, hjust = 1, color="black"), axis.line = element_line(colour = "black"), panel.background = element_rect(fill = "transparent", colour = NA), legend.position = "none")
