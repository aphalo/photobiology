library(lubridate)
# today() does not behave as I expected
tz(now(tzone = "Europe/Helsinki"))
tz(today(tzone = "Europe/Helsinki")) # UTC!!
today("America/New_York") == today("Asia/Tokyo") # FALSE with no warning.
tz(today("America/New_York")) # UTC!
tz(today("Asia/Tokyo")) # UTC!
ymd("2017-06-19", tz = "America/New_York") == ymd("2017-06-20", tz = "Asia/Tokyo") # FALSE with warning.
tz(ymd("2017-06-19", tz = "America/New_York")) # OK
tz(ymd("2017-06-20", tz = "Asia/Tokyo")) # OK

# nothing unusual here
tz(as.POSIXlt(now(), tz = "Europe/Helsinki"))
tz(as.POSIXlt(now(tzone = "Europe/Helsinki")))
tz(as.POSIXlt(now(tzone = "Europe/Helsinki"), tz = "Europe/Helsinki"))
# but the unexpected seems to be the behaviour of as.Date.POSIXlt() as all returned values have tz set to "UTC"
tz(as.Date(as.POSIXlt(now())))  # UTC!!
tz(as.Date(as.POSIXlt(now(tzone = "Europe/Helsinki"), tz = "Europe/Helsinki"), tz = "Europe/Helsinki"))  # UTC!!
tz(as.Date(as.POSIXlt(now(), tz = "Europe/Helsinki"), tz = "Europe/Helsinki"))  # UTC!!
tz(as.Date(as.POSIXlt(now()), tz = "Europe/Helsinki"))  # UTC!!
