
hist.marg.data <- scan('temp.hist-eqtl-marg.txt')
hist.dir.data <- scan('temp.hist-eqtl-dir.txt')


.data <- hist.dir.data
n <- length(.data)

.x <- -log10((1:n)/(n+1))
.y <- sort(.data, decreasing = TRUE)

plot(.x, .y)
abline(a=0,b=1)
