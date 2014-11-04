test.df <- data.frame(a=1:10, b=rep(1,10))
class(test.df) <- c("deriv", class(test.df))
class(test.df)
test0.df <- test.df
class(test0.df)
test1.df <- test.df[ test.df$a>=2 & test.df$a<=9 , ]
class(test1.df)
test2.df <- test.df[ 2:9, ]
class(test2.df)


library(data.table)
test.dt <- data.table(a=1:10, b=rep(1,10))
setattr(test.dt, "class", c("deriv", class(test.dt)))
class(test.dt)
test0.dt <- copy(test.dt)
class(test0.dt)
test1.dt <- test.dt[ a>=2 & a<=9 ]
class(test1.dt)
test2.dt <- test.dt[ a %between% c(2,9)]
class(test2.dt)
test3.dt <- test.dt[ 2:9 ]
class(test3.dt)

test.dt[ , b := NULL]
(test.dt[ , b :=b + 1])
test1.dt[ , b := NULL]
test2.dt[ , b := NULL]
test3.dt[ , b := NULL]

