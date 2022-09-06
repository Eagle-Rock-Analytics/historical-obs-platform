## https://drive.google.com/file/d/1tV0naJXM89RcJXoThdmdMi4WCD5QFBO_/view
## Steps to run this script:
## Place gzipped input data downloaded form main archive into a folder:
## d:/solar data/2014 solar data/ .
## Have folder ready:
## d:/solar data parsed/
## Issue these commands in R:
## filedate <- '20150101'
## parsecwoptxt(filedate)

parsecwoptxt <- function(filedate) {

  thisyear<-substr(filedate,1,4)
  thismonth<-substr(filedate,5,6)
  thisday<-substr(filedate,7,8)
  thisdate<-paste(thisyear,thismonth,thisday,sep="-")
  thisdateasdate<-as.Date(thisdate)

  yesterdayasdate<-thisdateasdate-1
  yesterdayyear<-substr(as.character(yesterdayasdate),1,4)
  yesterdaymonth<-substr(as.character(yesterdayasdate),6,7)
  yesterdayday<-substr(as.character(yesterdayasdate),9,10)
  yesterdaydate<-paste(yesterdayyear,yesterdaymonth,yesterdayday,sep="-")

  filename<-paste("d:/a solar data/2014 solar data/L",filedate,".txt.gz",sep="")
  mydataout <- paste("d:/solar data parsed/LP",filedate,"P.txt",sep="")

  unlista<-readLines(filename)
  report<-unlista

  stationpat<-"([a-zA-Z0-9-]{1,})(>)"
  statns<- str_match(unlista,stationpat)
  stationname<- statns[,2]

  timepat<-"(>)([0-9]{1,6})([a-z])"
  dtimes<- str_match(unlista,timepat)

  datatimes<-dtimes [,3]

  recordday<-substr(datatimes,1,2)
  rightday<-recordday==thisday
  rightyesterday<-recordday==yesterdayday
  wrongday<-!(recordday==thisday|recordday==yesterdayday)

  dateflag<-integer()
  dateflag[rightday]<-1
  dateflag[rightyesterday]<-2
  dateflag[wrongday]<-3

  recordhour<-substr(datatimes,3,4)
  recordmin<-substr(datatimes,5,6)
  recordtime<-paste(recordhour,recordmin,sep=":")
  recordstamp<-paste(thisdate,recordtime,sep=" ")

  z <- strptime(recordstamp, "%Y-%m-%d %H:%M",tz='UTC')
  z[dateflag!=1]<-NA

  # Next latitude and longitude.
  latpat<-"([a-z -]{1,})([0-9]{1,4}.[0-9]{1,2})([N n S s]{1,1})"
  latparts<-str_match(unlista,latpat)
  lat<-latparts[,3]
  latd<-as.numeric(substr(lat,1,2))
  latmin<-as.numeric(substr(lat,3,7))
  latitude<-(latd+(latmin/60))
  latsign<-latparts[,4]

  lonpat<-"(/)([0-9]{1,5}.[0-9]{1,6})([E e W w]{1,1})"
  lonparts<-str_match(unlista,lonpat)
  lon<-lonparts[,3]
  lond<-as.numeric(substr(lon,1,3))
  lonmin<-as.numeric(substr(lon,4,8))
  lonsign<-lonparts[,4]
  longitude<-(lond+(lonmin/60))

  # Next the data.
  windpat<-"(_)([0-9]{1,3})/([0-9]{1,3})"
  windparts<-str_match(unlista,windpat)
  winddir<-windparts[,3]
  windknots<-windparts[,4]

  gustpat<-"(g)([0-9]{1,3})([a-z]){1,1}"
  gparts<-str_match(unlista,gustpat)
  gust<-gparts[,3]

  tpat<-"(t)([- 0-9]{1,3})([a-z A-Z]){1,1}"
  tparts<-str_match(unlista,tpat)
  temp<-tparts[,3]

  rpat<-"(r)([0-9]{1,3})([a-z A-Z]){1,1}"
  rparts<-str_match(unlista,rpat)
  rainfallhour<-rparts[,3]

  ppat<-"(p)([0-9]{1,3})([a-z A-Z]){1,1}"
  pparts<-str_match(unlista,ppat)
  rainfall24h<-pparts[,3]

  bigppat<-"(P)([0-9]{1,3})"
  bigpparts<-str_match(unlista,bigppat)
  rainfalltoday<-bigpparts[,3]

  hpat<-"(h)([0-9]{1,2})"
  hparts<-str_match(unlista,hpat)
  relativehumidity<-hparts[,3]

  bpat<-"(b)([0-9]{1,5})"
  bparts<-str_match(unlista,bpat)
  baropressure<-bparts[,3]

  lpat<-"([L l])([0-9]{3,4})"
  lparts<-str_match(unlista,lpat)
  lrec<-lparts[,1]
  lcharerr<-nchar(lparts[,1])>4
  lrecord<-as.integer(lparts[,3])
  lstyle<-(lparts[,2]=="l")
  laddend<-lstyle*1000
  lfin<-lrecord+laddend

  techbreak<-str_locate(unlista,lpat)
  totalchar<-nchar(unlista)
  tsuffix<-totalchar-techbreak[,2]
  tech<-substrRight(unlista,tsuffix)

  outtable<-data.frame(report,stationname,thisdate,
  datatimes,dateflag,z,latitude,latsign,longitude,
  lonsign,winddir,windknots,gust,temp, rainfallhour,
  rainfall24h,rainfalltoday,relativehumidity,baropressure,
  lrec,lfin,lcharerr,tech)

  write.table(outtable,mydataout,sep="\t",row.names=FALSE)
}


# The output variables are as follows:
# report = the packet received at CWOP
# stationname = reported station name
# thisdate = the date of the data archive, expressed as a date
# datatimes = the timestamp within the packet
# z = the imputed time of the data expressed as a date
# latitude = absolute value of latitude
# latsign = N or S
# longitude = absolute value of the longitude
# lonsing = E or W
# windknots = wind in knots
# winddir = wind direction in degrees clockwise from north.
# gust = peak wind speed in mph in the last 5 minutes.
# temp = temperature (in degrees Fahrenheit (below zero as -01 to -99)
# rainfallhour = rainfall (in hundredths of an inch) in the last hour.
# rainfall24h = rainfall (in hundredths of an inch) in the last 24 hours.
# rainfalltoday = rainfall (in hundredths of an inch) since midnight.
# relativelhumidity = humidity (in %. 00 = 100%)
# baropressure = p in tenths of millibars/tenths of hPascal
# lrec = the luminosity record pre-parsing
# lfin = calculation of luminosity from that record
# lcharerr = a flag that is TRUE for protocols Lxxxx or lxxxx
# tech = the suffix that is supposed to describe tecnology.
# Note that the APRS protocol is Lxxx or lxxx: L = luminosity (in watts per
# square meter) 999 and below. while l, that is lower-case letter "L" =
# luminosity (in watts per square meter) 1000 and above.
