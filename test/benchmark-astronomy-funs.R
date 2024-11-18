# This are just benchmarks
# We need to check also if values are correct
library(microbenchmark)
library(svglite)
library(tibble)
library(dplyr)
library(ggplot2)
library(ggpmisc)
library(lubridate)
library(photobiology) # sun_angles() and day_night()
library(solartime) # computeSunPosition(), computeSunriseHour(), computeSunsetHour()
library(suncalc) # getSunlightPosition(), getSunlightTimes()
library(fishmethods) # astrocalc4r()

## From Gary Nelson, maintainer of 'fishmethods'

astrocalc4r_Version2.3_for_fishmethods=function (day, month, year, hour, timezone, lat, lon, withinput = FALSE,
                                                 seaorland = "maritime", acknowledgment = FALSE)
{
  if (acknowledgment) {
    cat("\n", "---------------------------------------------------------")
    cat("\n", "                AstroCalcPureR Version 2.3")
    cat("\n", "Documentation: Jacobson L, Seaver A, Tang J. 2011. AstroCalc4R:")
    cat("\n", "software to calculate solar zenith angle; time at sunrise,")
    cat("\n", "local noon and sunset; and photosynthetically available")
    cat("\n", "radiation based on date, time and location. US Dept Commer,")
    cat("\n", "Northeast Fish Sci Cent Ref Doc. 11-14; 10 p. Available from:")
    cat("\n", "National Marine Fisheries Service, 166 Water Street, ")
    cat("\n", "Woods Hole, MA 02543-1026, or online at")
    cat("\n", "http://nefsc.noaa.gov/publications/")
    cat("\n \n", "Available in fishmethods library.  Contact the fishmethods")
    cat("\n" , "administrator or Larry Jacobson (NOAA, National Marine")
    cat("\n", "Fisheries Service-retired) at larryjacobson6@gmail.com")
    cat("\n", "for assitance.")
    cat("\n\n", "Useage:")
    cat("\n", "    AstroCalcPureR(day,month,year,hour,timezone,")
    cat("\n", "                   lat,lon,withinput=F,")
    cat("\n", "                   seaorland='maritime',")
    cat("\n", "                   acknowledgment=TRUE)")
    cat("\n\n", "HINT: set acknowledgment=FALSE to avoid this message")
    cat("\n", "---------------------------------------------------------",
        "\n")
  }
  ###
  # Modification history:
  # Added traps to ensure that timezone and longitude have the same sign if
  # time is not UTC (timezone==0).  In particular, if time is not UTC
  # then both should be negative in western hemisphere or both positive in
  # the eastern hemisphere.  This is a bit tricky because
  # the boundaries for time zones are not based strictly on latitude.  The
  # solution here is to require UTC (universal stardard time inputs) input data
  # (argument timezone=0) if  there is a conflict.  In other words, require
  # the user to convert input time to UTC and use timezone=0 if necessary.
  #
  # Some minor changes to error messages.
  #
  # Version number in optional acknoledgement changed from to 2.2 to 2.3.
  #                        - Larry Jacobson, January 13, 2019
  #
  # timezone=0 if a sign conflict
  options(digits = 9)
  deg2rad <- pi/180
  null.c <- function(x) return(sum(is.null(x)))
  if (sum(null.c(day), null.c(month), null.c(year), null.c(hour),
          null.c(timezone), null.c(lat), null.c(lon)) > 0)
    stop("\n Null or missing required data vector for day, month, year, timezone, lat or lon \n")
  if ((length(day) != length(month)) | (length(month) != length(year)) |
      (length(year) != length(hour)) | (length(hour) != length(timezone)) |
      (length(timezone) != length(lat)) | (length(lat) !=
                                           length(lon)))
    stop("\n Input vectors are not the same length \n")
  times <- length(day)
  na.c <- function(x) return(sum(is.na(x)))
  if (sum(na.c(day), na.c(month), na.c(year), na.c(hour),
          na.c(timezone), na.c(lat), na.c(lon)) > 0)
    stop("\n NA values in input data \n")
  logic1 <- year < 0
  if (sum(logic1) > 0)
    stop(cat("\n Error in year at rows:", (1:times)[logic1],
             " \n\n"))
  is.leap <- function(x) return((((x%%4 == 0) & (x%%100 !=
                                                   0))) | (x%%400 == 0))
  date.list <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30,
                 31)
  logic1 <- abs(month - 6) > 6
  if (sum(logic1) > 0)
    stop(cat("\n Error in month at rows:", (1:times)[logic1],
             " \n\n"))
  logic1 <- day > (date.list[month] + is.leap(year) * (month ==
                                                         2))
  logic2 <- day <= 0
  if ((sum(logic1) > 0) | (sum(logic2) > 0))
    stop(cat("\n Incorrect month-day-year combination at rows: ",
             (1:times)[logic1 | logic2], " \n\n"))
  logic1 <- abs(hour - 12) > 12
  if (sum(logic1) > 0)
    stop(cat("\n Error in hour at rows:", (1:times)[logic1],
             " \n\n"))
  logic1 <- abs(timezone) > 12
  if (sum(logic1) > 0)
    stop(cat("\n Error in time zone at rows:",
             (1:times)[logic1], " \n\n"))
  logic1 <- abs(lon) > 180
  if (sum(logic1) > 0)
    stop(cat("\n Error in longitude at rows:",
             (1:times)[logic1], " \n\n"))
  logic1 <- abs(lat) > 90
  if (sum(logic1) > 0)
    stop(cat("\n Error in latitude at rows:",
             (1:times)[logic1], " \n\n"))
  logic1 <- sign(lon) == sign(timezone)
  logic2 <- timezone == 0
  logic3 <- !(logic1 | logic2)
  #  print(c(logic1=logic1,logic2=logic2,logic3))
  if(sum(logic3) !=0) stop(cat("\n \n Arguments longitude and timezone must have the same sign if input time is",
                               "\n not UTC (timezone != 0).  In particular, if timezone !=0, both lon and timezone must",
                               "\n be negative for locations in western hemisphere and positive for locations in the",
                               "\n eastern hemisphere.  Check and fix input data if warranted. If data are correct",
                               "\n then convert input time (argument hour) to UTC and use timezone=zero.",
                               "\n This problem  occurs ",sum(logic3)," times at rows: ",(1:times)[logic3],"\n\n"))

  JulianDay <- function(xday, xmonth, xyear) {
    mm <- xmonth
    xmonth[mm <= 2] <- xmonth[mm <= 2] + 12
    xyear[mm <= 2] <- xyear[mm <= 2] - 1
    xa <- floor(xyear/100)
    xb <- 2 - xa + floor(xa/4)
    jd <- floor(365.25 * (xyear + 4716)) + floor(30.6001 *
                                                   (xmonth + 1)) + xday + xb - 1524.5
    return(jd)
  }
  daymonth <- function(mth, yr) {
    day[is.leap(yr)] <- c(31, 29, 31, 30, 31, 30, 31, 31,
                          30, 31, 30, 31)[mth[is.leap(yr)]]
    day[!is.leap(yr)] <- c(31, 28, 31, 30, 31, 30, 31, 31,
                           30, 31, 30, 31)[mth[!is.leap(yr)]]
    return(day)
  }
  parcalc <- function(zenith, setting = seaorland) {
    I0 <- 531.2
    V <- 23
    uv <- 1.4
    u0 <- 0.34
    r <- 0.05
    d <- 1
    if (!setting %in% c("maritime", "continental"))
      stop("setting value is neither 'maritime' nor 'continental'!")
    if (setting == "maritime") {
      a <- 0.068
      b <- 0.379
      a1 <- 0.117
      b1 <- 0.493
      av <- 0.002
      bv <- 0.87
      a0 <- 0.052
      b0 <- 0.99
    }
    else if (setting == "continental") {
      a <- 0.078
      b <- 0.882
      a1 <- 0.123
      b1 <- 0.594
      av <- 0.002
      bv <- 0.87
      a0 <- 0.052
      b0 <- 0.99
    }
    zrad <- zenith * deg2rad
    x1 <- uv/cos(zrad)
    xx <- exp(-av * x1^bv)
    x2 <- u0/cos(zrad)
    xxx <- exp(-a0 * x2^b0)
    xa <- a + b/V
    xb <- d - r * (a1 + b1/V)
    par <- I0 * cos(zrad) * exp(-xa/cos(zrad))/xb * xx *
      xxx
    par[zenith > 89.9999] <- 0
    return(par)
  }
  output <- as.data.frame(matrix(nrow = 0, ncol = 9))
  names(output) <- c("noon", "sunrise", "sunset", "azimuth",
                     "zenith", "eqtime", "declin", "daylight", "PAR")
  hourtemp <- hour - timezone
  hour <- ifelse(hourtemp > 24, hourtemp - 24, hourtemp)
  change_day <- !(hour == hourtemp)
  dm <- daymonth(month, year)
  daytemp <- day
  daytemp[change_day] <- ifelse((day[change_day] < dm[change_day]),
                                day[change_day] + 1, 1)
  change_month <- abs(day - daytemp) > 1
  monthtemp <- month
  monthtemp[change_month] <- ifelse(month[change_month] <
                                      12, month[change_month] + 1, 1)
  change_year <- abs(month - monthtemp) > 1
  yeartemp <- year
  yeartemp[change_year] <- year[change_year] + 1
  xy <- yeartemp
  xm <- monthtemp
  xd <- daytemp + hourtemp/24
  jd <- JulianDay(xd, xm, xy) * 100/100
  jc <- (jd - 2451545)/36525
  xx <- 280.46646 + 36000.76983 * jc + 0.0003032 * jc^2
  gmls <- xx%%360
  xx <- 357.52911 + 35999.05029 * jc - 0.0001537 * jc^2
  gmas <- xx%%360
  eeo <- 0.016708634 - 4.2037e-05 * jc - 1.267e-07 * jc^2
  scx <- (1.914602 - 0.004817 * jc - 1.4e-05 * jc^2) * sin(gmas *
                                                             deg2rad) + (0.019993 - 0.000101 * jc) * sin(2 * gmas *
                                                                                                           deg2rad) + 0.000289 * sin(3 * gmas * deg2rad)
  Stl <- gmls + scx
  Sta <- gmas + scx
  srv <- 1.000001018 * (1 - eeo^2)/(1 + eeo * cos(Sta * deg2rad))
  omega <- 125.04 - 1934.136 * jc
  lambda <- Stl - 0.00569 - 0.00478 * sin(omega * deg2rad)
  epsilon <- (23 + 26/60 + 21.448/60^2) - (46.815/60^2) *
    jc - (0.00059/60^2) * jc^2 + (0.001813/60^2) * jc^3
  oblx <- 0.00256 * cos(omega * deg2rad)
  epsilon <- epsilon + oblx
  alpha <- atan2(cos(epsilon * deg2rad) * sin(lambda * deg2rad),
                 cos(lambda * deg2rad))/deg2rad
  declin <- asin(sin(epsilon * deg2rad) * sin(lambda * deg2rad))/deg2rad
  y <- tan(epsilon * deg2rad/2)^2
  eqtime <- (y * sin(2 * gmls * deg2rad) - 2 * eeo * sin(gmas *
                                                           deg2rad) + 4 * eeo * y * sin(gmas * deg2rad) * cos(2 *
                                                                                                                gmls * deg2rad) - y^2 * sin(4 * gmls * deg2rad)/2 -
               5/4 * eeo^2 * sin(2 * gmas * deg2rad))/deg2rad * 4
  h0 <- -0.8333 * deg2rad
  phi <- lat * deg2rad
  hangle <- acos((sin(h0) - sin(declin * deg2rad) * sin(phi))/cos(declin *
                                                                    deg2rad)/cos(phi))/deg2rad
  noon <- (720 - 4 * lon + timezone * 60 - eqtime)/1440
  sunrise <- (noon * 1440 - hangle * 4)/1440 * 24
  sunset <- (noon * 1440 + hangle * 4)/1440 * 24
  noon <- noon * 24
  daylight <- hangle * 8
  tst <- (hourtemp * 60 + eqtime + 4 * lon)%%1440
  tsa <- ifelse(tst < 0, tst/4 + 180, tst/4 - 180)
  zenith <- 90 - asin(sin(lat * deg2rad) * sin(declin * deg2rad) +
                        cos(lat * deg2rad) * cos(declin * deg2rad) * cos(tsa *
                                                                           deg2rad))/deg2rad
  azimuth <- acos((sin(lat * deg2rad) * sin((90 - zenith) *
                                              deg2rad) - sin(declin * deg2rad))/cos(lat * deg2rad)/cos((90 -
                                                                                                          zenith) * deg2rad))/deg2rad + 180
  azimuth <- ifelse(tsa > 0, azimuth%%360, 360 - azimuth%%360)
  daylight <- daylight/60
  PAR <- parcalc(zenith)
  if (any(is.nan(sunrise))) {
    message(paste("Warning: Polar day/night (daylength 0 or 24 hrs) at record(s):",
                  (1:times)[is.nan(sunrise)], "\n Check input data (i.e. latitude)?"))
    daylight <- ifelse(PAR > 0, 24, 0)
  }
  output <- rbind(output, data.frame(noon = noon, sunrise = sunrise,
                                     sunset = sunset, azimuth = azimuth, zenith = zenith,
                                     eqtime = eqtime, declin = declin, daylight = daylight,
                                     PAR = PAR))
  if (withinput)
    return(cbind(data.frame(tzone = timezone, day = day,
                            month = month, year = year, hhour = hour, xlat = lat,
                            xlon = lon), output))
  else return(output)
}

## vectorised time points

num.years <- 1
times.per.min <- as.POSIXct(seq(from = today(tzone = "Europe/Helsinki"), to = today(tzone = "Europe/Helsinki") + years(num.years),
                 length.out =  365 * 24 * 60 * (num.years)))
times.per.hour <- as.POSIXct(seq(from = today(tzone = "Europe/Helsinki"), to = today(tzone = "Europe/Helsinki") + years(num.years),
                             length.out =  365 * 24 * (num.years)))
times.per.day <- as.POSIXct(seq(from = today(tzone = "Europe/Helsinki"), to = today(tzone = "Europe/Helsinki") + years(num.years),
                             length.out =  365 * (num.years)))
times.per.month <- as.POSIXct(seq(from = today(tzone = "Europe/Helsinki"), to = today(tzone = "Europe/Helsinki") + years(num.years),
                            length.out =  12 * (num.years)))
times.per.year <- as.POSIXct(seq(from = today(tzone = "Europe/Helsinki"), to = today(tzone = "Europe/Helsinki") + years(num.years),
                             length.out = num.years))
length(times.per.year)
times.per.year

times <- list(per.min = times.per.min, per.hour = times.per.hour, per.day = times.per.day, per.month = times.per.month, per.year = times.per.year)

benchmark.results <- list()
idx <- 0L
for (t in times) {
  idx <- idx + 1L
  num.evals <- as.integer(max(7, min(21, length(times.per.hour) / length(t))))
  result <-
    microbenchmark(sun_angles(t, geocode = data.frame(lat = 60, lon = 0)),
                   computeSunPosition(timestamp = t, 60, 0),
                   times = num.evals, unit = "s",
                   setup = gc())
  result$case <- names(times)[idx]
  result$case.length <- length(t)
  benchmark.results <- c(benchmark.results, list(as.data.frame(result)))
  print(result)
    cat(".")
}
result.vec.times.tb <- bind_rows(benchmark.results)
bind_rows(benchmark.results) %>%
  mutate(fun = ifelse(grepl("sun_angles", expr), "photobiology::sun_angles", "solartime::computeSunPosition")) %>%
  group_by(case, fun) %>%
  summarise(case.length = median(case.length),
            median.time = median(time) * 1e-9, # seconds
            median.time.per.time.point = median.time / case.length) -> summary.vec.times.tb
summary.vec.times.tb

summary.vec.times.tb %>%
  ggplot(aes(case.length, median.time * 1000, colour = fun)) +
  geom_point() +
  geom_line() +
  scale_x_log10(name = "Length of vector of times points") +
  scale_y_log10(name = "Total execution time (ms)") +
  theme_bw() + theme(legend.position = "top") -> total.times.ggp

svglite("./test/total-times-ggp.svg", width = 6, height = 4)
print(total.times.ggp)
dev.off()

summary.vec.times.tb %>%
  ggplot(aes(case.length, median.time.per.time.point * 1e6, colour = fun)) +
  geom_point() +
  geom_line() +
  expand_limits(y = 0) +
  scale_x_log10(name = "Length of vector of time points") +
  scale_y_continuous(name = expression("Execution time per time point"~~(mu*s))) +
  theme_bw() + theme(legend.position = "top") -> per.timepoint.ggp

svglite("./test/per-timepoint-ggp.svg", width = 6, height = 4)
print(per.timepoint.ggp)
dev.off()

## Vectorised latitudes

latitudes <- list(thousands = seq(0, 87, length.out = 2^12),
                  hundreds = seq(0, 87, length.out = 2^9),
                  tens = seq(0, 87, length.out = 2^6),
                  ones = seq(0, 87, length.out = 2^3))

# every hour for a week
t <- as.POSIXct(seq(from = today(tzone = "Europe/Helsinki"), to = today(tzone = "Europe/Helsinki") + days(1),
                    length.out =  24))

benchmark.results <- list()
idx <- 0L
for (l in latitudes) {
  idx <- idx + 1L
  num.evals <- as.integer(max(11, min(99, 5 * length(latitudes[[1]]) / length(l))))
  localities <- data.frame(lat = l, lon = 0)
  result <-
    microbenchmark(sun_angles(t, geocode = localities),
                   for (i in 1:nrow(localities)) {
                     with(localities[i, ],
                          computeSunPosition(timestamp = t, lat, lon)
                     )
                   },
                   times = num.evals, unit = "s",
                   setup = gc())
  result$case <- names(latitudes)[idx]
  result$case.length <- length(l) * length(t)
  result$num.localities <- length(l)
  benchmark.results <- c(benchmark.results, list(as.data.frame(result)))
  print(result)
  cat(".")
}
result.vec.latitudes.tb <- bind_rows(benchmark.results)
result.vec.latitudes.tb %>%
  mutate(fun = ifelse(grepl("sun_angles", expr), "photobiology::sun_angles", "solartime::computeSunPosition")) %>%
  group_by(case, fun) %>%
  summarise(case.length = median(case.length),
            num.localities = median(num.localities),
            median.time = median(time) * 1e-9, # seconds
            median.time.per.time.point = median.time / case.length) -> summary.vec.latitudes.tb
summary.vec.latitudes.tb

summary.vec.latitudes.tb %>%
  ggplot(aes(num.localities, median.time * 1000, colour = fun)) +
  geom_point() +
  geom_line() +
  scale_x_log10(name = "Length of vector of latitudes") +
  scale_y_log10(name = "Total execution time (ms)") +
  theme_bw() + theme(legend.position = "top") -> total.times.latitudes.ggp

svglite("./test/total-times-latitudes-ggp.svg", width = 6, height = 4)
print(total.times.latitudes.ggp)
dev.off()

summary.vec.latitudes.tb %>%
  ggplot(aes(num.localities, median.time.per.time.point * 1e6, colour = fun)) +
  geom_point() +
  geom_line() +
  scale_x_log10(name = "Length of vector of latitudes") +
  scale_y_continuous(name = expression("Execution time per time point"~~(mu*s))) +
  expand_limits(y = 0) +
  theme_bw() + theme(legend.position = "top") -> per.timepoint.latitudes.ggp

svglite("./test/per-timepoint-latitudes-ggp.svg", width = 6, height = 4)
print(per.timepoint.latitudes.ggp)
dev.off()

## profile

t <- as.POSIXct(seq(from = today(tzone = "Europe/Helsinki"), to = today(tzone = "Europe/Helsinki") + days(1),
                    length.out = 3))
localities <- data.frame(lat = latitudes$thousands, lon = 0)
profvis::profvis(sun_angles(t, geocode = localities))

## test values computed for extreme dates

date_times <-  as.POSIXct((ymd_hm("0001-01-01 00:00") +
                             years(c(0, 100, 1000, 1800, 1900, 2000, 2020, 2100, 3000) - 1)) +
                            days(rep((0:11) * 30, 9)) +
                            hours(rep(c(0,6,12,18), 9 * 12)))
localities <- data.frame(lat = latitudes$tens, lon = 0)[-1, ] # computeSunPosition fails at the equator

length(date_times) * nrow(localities)

## this need to be expanded to multiple times and locations
time <- date_times[7] + months(2)
astrocalc4r(day = day(time), month = month(time), year = year(time), hour = hour(time) + 0.06,
             timezone = 0, lat = localities[50, "lat"], localities[50, "lon"],
            withinput = TRUE, seaorland = "continental")

day_night(time + minutes(3) + seconds(40), geocode = localities[50, ])

sun_angles(time + minutes(3) + seconds(40), geocode = localities[50, ])

getSunlightPosition(time + minutes(3) + seconds(40),
                    lat = localities[50, "lat"], lon = localities[50, "lon"]) %>%
  mutate(altitude = altitude * 180 / pi,
         azimuth = (pi + azimuth) * 180 / pi %% 360)

## photobiology

sa_results.tb <- cbind(sun_angles(date_times, geocode = localities),
                           dplyr::select(day_night(date_times, geocode = localities),
                                  -day, -tz, -longitude, -latitude, -address))
names(sa_results.tb) <- paste("sa.", names(sa_results.tb), sep = "")

nrow(sa_results.tb)

###

csp_results.lst <- list(length(localities))
sum_nrows <- 0
for (i in 1:nrow(localities)) {
  with(localities[i, ],
       computeSunPosition(date_times, lat, lon)) %>%
    as.data.frame() %>%
    mutate(declination = declination * 180 / pi,
           elevation = elevation * 180 / pi,
           azimuth = azimuth * 180 / pi) -> tmp
  sum_nrows <- sum_nrows + nrow(tmp)
  csp_results.lst[[i]] <- tmp
}
print(sum_nrows)

csp_results.tb <- dplyr::bind_rows(csp_results.lst)
names(csp_results.tb) <- paste("csp.", names(csp_results.tb), sep = "")

##

sc_results.lst <- list(length(localities))
sum_nrows <- 0
for (i in 1:nrow(localities)) {
  with(localities[i, ],
       getSunlightPosition(date_times,
                           lat = lat, lon = lon)) %>%
         mutate(altitude = altitude * 180 / pi,
                azimuth = (pi + azimuth) * 180 / pi %% 360) -> tmp
  sum_nrows <- sum_nrows + nrow(tmp)
  sc_results.lst[[i]] <- tmp
}
print(sum_nrows)

sc_results.tb <- dplyr::bind_rows(sc_results.lst)
names(sc_results.tb) <- paste("sc.", names(sc_results.tb), sep = "")

##

ac4r_results.lst <- list(length(localities) * length(date_times))
sum_nrows <- 0
for (i in 1:nrow(localities)) {
  for (j in 1:length(date_times)) {
    astrocalc4r(day = day(date_times[j]), month = month(date_times[j]), year = year(date_times[j]), hour = hour(date_times[j]),
                timezone = 0, lat = localities[["lat"]][1], lon = localities[["lon"]][1],
                withinput = TRUE, seaorland = "continental") -> tmp
  sum_nrows <- sum_nrows + nrow(tmp)
  ac4r_results.lst[[(i - 1) * length(date_times) + j]] <- tmp
  }
}
print(sum_nrows)

ac4r_results.tb <- dplyr::bind_rows(ac4r_results.lst)
names(ac4r_results.tb) <- paste("acr.", names(ac4r_results.tb), sep = "")

##

all_results.tb <- bind_cols(sa_results.tb, csp_results.tb, sc_results.tb, ac4r_results.tb) %>%
  mutate(sa_year = year(sa.time),
         sa_hour = hour(sa.time),
         sa_hour_solar = round(sa.solartime),
         sa_hour_solar = ifelse(sa_hour_solar == 24, 0, sa_hour_solar),csp_delta_elevation = csp.elevation - sa.elevation,
         csp_delta_azimuth = csp.azimuth - sa.azimuth,
         # csp_delta_azimuth = ifelse(csp_delta_azimuth < -180, NA, csp_delta_azimuth),
         # csp_delta_azimuth = ifelse(csp_delta_azimuth >  180, NA, csp_delta_azimuth),
         csp_delta_elevation = csp.elevation - sa.elevation,
         csp.solartime = csp.hour,
         csp_delta_solar_time = ifelse(csp.solartime < 0, 24 - csp.solartime, csp.solartime) - sa.solartime,

         sc_delta_azimuth = sc.azimuth - acr.azimuth,
         sc_delta_elevation = sc.altitude - sa.elevation,

         acr_delta_azimuth = acr.azimuth - sa.azimuth,
         acr.elevation = -(acr.zenith - 90),
         acr_delta_elevation = acr.elevation - sa.elevation,
         acr_delta_declination = acr.declin - sa.declination,
         acr_delta_daylength = acr.daylight - sa.daylength,
         acr_delta_sunset = acr.sunset - sa.sunset,
         acr_delta_sunrise = acr.sunrise - sa.sunrise,
         acr_delta_noon = acr.noon - sa.noon
  )

### Azimuth

all_results.tb %>%
  filter(sa_year != 3000 & sa_year >= 1000) %>%
  ggplot(aes(sa.latitude, acr_delta_azimuth, colour = factor(sa_hour))) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "red") +
  geom_point(alpha = 1/3) +
  stat_quadrant_counts(size = 3) +
  expand_limits(y = c(-2, 2)) +
  facet_wrap(~sa_year, labeller = label_both, ncol = 2) +
  labs(y = "Difference in solar azimuth (degrees)",
       x = "Latitude (degrees)",
       colour = "Time of day\nUTC (h)") +
  theme_bw() + theme(legend.position = "top") +
  ggtitle("Azimuth: astrocalc4r() vs. sun_angles()") -> acr4r_azimuth_error.ggp

svglite("./test/acr4r-azimuth-error-ggp.svg", width = 6, height = 10)
print(acr4r_azimuth_error.ggp)
dev.off()

all_results.tb %>%
  filter(sa_year != 3000 & sa_year >= 1000) %>%
  ggplot(aes(sa.latitude, csp_delta_azimuth, colour = factor(sa_hour))) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "red") +
  geom_point(alpha = 1/3) +
  stat_quadrant_counts(size = 3) +
  expand_limits(y = c(-2, 2)) +
  facet_wrap(~sa_year, labeller = label_both, ncol = 2) +
  labs(y = "Difference in solar azimuth (degrees)",
       x = "Latitude (degrees)",
       colour = "Time of day\nUTC (h)") +
  theme_bw() + theme(legend.position = "top") +
  ggtitle("Azimuth: computeSunPosition() vs. sun_angles()") -> csp_azimuth_error.ggp

svglite("./test/csp-azimuth-error-ggp.svg", width = 6, height = 10)
print(csp_azimuth_error.ggp)
dev.off()


all_results.tb %>%
  filter(sa_year != 3000 & sa_year >= 1000) %>%
  ggplot(aes(sa.latitude, sc_delta_azimuth, colour = factor(sa_hour))) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "red") +
  geom_point(alpha = 1/3) +
  stat_quadrant_counts(size = 3) +
  expand_limits(y = c(-2, 2)) +
  facet_wrap(~sa_year, labeller = label_both, ncol = 2) +
  labs(y = "Difference in solar azimuth (degrees)",
       x = "Latitude (degrees)",
       colour = "Time of day\nUTC (h)") +
  theme_bw() + theme(legend.position = "top") +
  ggtitle("Azimuth: getSunlightPosition() vs. sun_angles()") -> sc_azimuth_error.ggp

svglite("./test/sc-azimuth-error-ggp.svg", width = 6, height = 10)
print(sc_azimuth_error.ggp)
dev.off()

### Elevation

all_results.tb %>%
  filter(sa_year != 3000 & sa_year >= 1000) %>%
  ggplot(aes(sa.latitude, acr_delta_elevation, colour = factor(sa_hour))) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "red") +
  geom_point(alpha = 1/3) +
  stat_quadrant_counts(size = 3) +
  expand_limits(y = c(-2, 2)) +
  facet_wrap(~sa_year, labeller = label_both, ncol = 2) +
  labs(y = "Difference in solar elevation (degrees)",
       x = "Latitude (degrees)",
       colour = "Time of day\nUTC (h)") +
  theme_bw() + theme(legend.position = "top") +
  ggtitle("Azimuth: astrocalc4r() vs. sun_angles()") -> acr4r_elevation_error.ggp

svglite("./test/acr4r-elevation-error-ggp.svg", width = 6, height = 10)
print(acr4r_elevation_error.ggp)
dev.off()

all_results.tb %>%
  filter(sa_year != 3000 & sa_year >= 1000) %>%
  ggplot(aes(sa.latitude, csp_delta_elevation, colour = factor(sa_hour))) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "red") +
  geom_point(alpha = 1/3) +
  stat_quadrant_counts(size = 3) +
  expand_limits(y = c(-2, 2)) +
  facet_wrap(~sa_year, labeller = label_both, ncol = 2) +
  labs(y = "Difference in solar elevation (degrees)",
       x = "Latitude (degrees)",
       colour = "Time of day\nUTC (h)") +
  theme_bw() + theme(legend.position = "top") +
  ggtitle("Azimuth: computeSunPosition() vs. sun_angles()") -> csp_elevation_error.ggp

svglite("./test/csp-elevation-error-ggp.svg", width = 6, height = 10)
print(csp_elevation_error.ggp)
dev.off()


all_results.tb %>%
  filter(sa_year != 3000 & sa_year >= 1000) %>%
  ggplot(aes(sa.latitude, sc_delta_elevation, colour = factor(sa_hour))) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "red") +
  geom_point(alpha = 1/3) +
  stat_quadrant_counts(size = 3) +
  expand_limits(y = c(-2, 2)) +
  facet_wrap(~sa_year, labeller = label_both, ncol = 2) +
  labs(y = "Difference in solar elevation (degrees)",
       x = "Latitude (degrees)",
       colour = "Time of day\nUTC (h)") +
  theme_bw() + theme(legend.position = "top") +
  ggtitle("Azimuth: getSunlightPosition() vs. sun_angles()") -> sc_elevation_error.ggp

svglite("./test/sc-elevation-error-ggp.svg", width = 6, height = 10)
print(sc_elevation_error.ggp)
dev.off()


#### Solar time

all_results.tb %>%
  filter(sa_year != 3000 & sa_year >= 1000) %>%
  ggplot(aes(sa.latitude, csp_delta_solar_time * 60, colour = factor(sa_hour))) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "red") +
  geom_point(alpha = 1/3) +
  stat_quadrant_counts(size = 3) +
  facet_wrap(~sa_year, labeller = label_both, ncol = 2) +
  labs(y = "Difference in solar time (min)",
       x = "Latitude (degrees)",
       colour = "Time of day\nUTC (h)") +
  theme_bw() + theme(legend.position = "top") +
  ggtitle("Azimuth: computeSunPosition() vs. sun_angles()") -> acr_solartime_error.ggp

svglite("./test/acr-solartime-error-ggp.svg", width = 6, height = 10)
print(acr_solartime_error.ggp)
dev.off()

### Sunrise

all_results.tb %>%
  filter(sa_year != 3000 & sa_year >= 1000) %>%
  ggplot(aes(sa.latitude, acr_delta_sunrise * 60, colour = factor(sa_hour))) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "red") +
  geom_point(alpha = 1/3) +
  stat_quadrant_counts(size = 3) +
  expand_limits(y = c(-2, 2)) +
  facet_wrap(~sa_year, labeller = label_both, ncol = 2) +
  labs(y = "Difference in sunrise time (min)",
       x = "Latitude (degrees)",
       colour = "Time of day\nUTC (h)") +
  theme_bw() + theme(legend.position = "top") +
  ggtitle("Azimuth: astrocalc4r() vs. sun_angles()") -> acr4r_sunrise_error.ggp

svglite("./test/acr4r-sunrise-error-ggp.svg", width = 6, height = 10)
print(acr4r_sunrise_error.ggp)
dev.off()

### Noon

all_results.tb %>%
  filter(sa_year != 3000 & sa_year >= 1000) %>%
  ggplot(aes(sa.latitude, acr_delta_noon * 60, colour = factor(sa_hour))) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "red") +
  geom_point(alpha = 1/3) +
  stat_quadrant_counts(size = 3) +
  expand_limits(y = c(-2, 2)) +
  facet_wrap(~sa_year, labeller = label_both, ncol = 2) +
  labs(y = "Difference in noon time (min)",
       x = "Latitude (degrees)",
       colour = "Time of day\nUTC (h)") +
  theme_bw() + theme(legend.position = "top") +
  ggtitle("Azimuth: astrocalc4r() vs. sun_angles()") -> acr4r_noon_error.ggp

svglite("./test/acr4r-noon-error-ggp.svg", width = 6, height = 10)
print(acr4r_noon_error.ggp)
dev.off()

### Sunset

all_results.tb %>%
  filter(sa_year != 3000 & sa_year >= 1000) %>%
  ggplot(aes(sa.latitude, acr_delta_sunset * 60, colour = factor(sa_hour))) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "red") +
  geom_point(alpha = 1/3) +
  stat_quadrant_counts(size = 3) +
  expand_limits(y = c(-2, 2)) +
  facet_wrap(~sa_year, labeller = label_both, ncol = 2) +
  labs(y = "Difference in sunset time (min)",
       x = "Latitude (degrees)",
       colour = "Time of day\nUTC (h)") +
  theme_bw() + theme(legend.position = "top") +
  ggtitle("Azimuth: astrocalc4r() vs. sun_angles()") -> acr4r_sunset_error.ggp

svglite("./test/acr4r-sunset-error-ggp.svg", width = 6, height = 10)
print(acr4r_sunset_error.ggp)
dev.off()

### Daylength

all_results.tb %>%
  filter(sa_year != 3000 & sa_year >= 1000) %>%
  ggplot(aes(sa.latitude, acr_delta_daylength * 60, colour = factor(sa_hour))) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "red") +
  geom_point(alpha = 1/3) +
  stat_quadrant_counts(size = 3) +
  expand_limits(y = c(-2, 2)) +
  facet_wrap(~sa_year, labeller = label_both, ncol = 2) +
  labs(y = "Difference in daylength time (min)",
       x = "Latitude (degrees)",
       colour = "Time of day\nUTC (h)") +
  theme_bw() + theme(legend.position = "top") +
  ggtitle("Azimuth: astrocalc4r() vs. sun_angles()") -> acr4r_daylength_error.ggp

svglite("./test/acr4r-daylength-error-ggp.svg", width = 6, height = 10)
print(acr4r_daylength_error.ggp)
dev.off()

### Actual values

all_results.tb %>%
  filter(sa_year != 3000 & sa_year >= 1000) %>%
  ggplot(aes(sa.latitude, colour = factor(sa_hour))) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "red") +
  geom_point(mapping = aes(y = sa.sunset), alpha = 1/3) +
  geom_point(mapping = aes(y = acr.sunset), alpha = 1/3, shape = "triangle") +
  expand_limits(y = c(-2, 2)) +
  facet_wrap(~sa_year, labeller = label_both, ncol = 2) +
  labs(y = "Sunset time (h)",
       x = "Latitude (degrees)",
       colour = "Time of day\nUTC (h)") +
  theme_bw() + theme(legend.position = "top") +
  ggtitle("Azimuth: astrocalc4r() vs. sun_angles()") -> sunset.ggp

svglite("./test/sunset-ggp.svg", width = 6, height = 10)
print(sunset.ggp)
dev.off()

## other functions

days <- seq(from = today(tzone = "Europe/Helsinki") - years(1000), to = today(tzone = "Europe/Helsinki"),
            length.out =  365 * 1000)
length(days)

microbenchmark(sunrise_times_photobiology.df <- day_night(days, geocode = data.frame(lat = 60, lon = 0)),
               times = 30, unit = "s")
sunrise_times_photobiology.df <- day_night(days, geocode = data.frame(lat = 60, lon = 0))

microbenchmark(sunrise_times_solartime.df <- computeSunriseHour(days, 60, 0),
               times = 30, unit = "s")
sunrise_times_solartime.df <- computeSunriseHour(days, 60, 0)

library(suncalc)

microbenchmark(sunrise_times_suncalc.df <- getSunlightTimes(days[1:(365*10)], 60, 0),
               times = 5, unit = "s")
sunrise_times_suncalc.df <- getSunlightTimes(days, 60, 0)

sun_rise.bench <-
  microbenchmark(getSunlightTimes(days[1:(365*50)], 60, 0),
                 computeSunriseHour(days[1:(365*50)], 60, 0),
                 day_night(days[1:(365*50)], geocode = data.frame(lat = 60, lon = 0), unit.out = "day"),
                 day_night(days[1:(365*50)], geocode = data.frame(lat = 60, lon = 0), unit.out = "datetime"),
                 times = 5, unit = "s")

autoplot(sun_rise.bench)

# test values

geocode <- data.frame(lat = 60.17, lon = 24.94, address = "Helsinki, Finland")
day_night(geocode = geocode, tz = "Europe/Helsinki", unit.out = "datetime")
day_night(geocode = geocode, tz = "Europe/Helsinki", unit.out = "hour")
local_t <- day_night(geocode = geocode, tz = "Europe/Helsinki", unit.out = "datetime")$sunrise
local_t
sol_t <- solar_time(local_t, geocode = geocode)
sol_t
computeSunriseHour(today(tzone = "Europe/Helsinki"), longDeg = geocode$lon, latDeg = geocode$lat)
