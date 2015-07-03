setMethod("plot", signature(x="pwm"), function(x, y="missing",...) {
  seqLogo(x)
})

setMethod("summary", signature(object="pwm"), function(object,...){
  cat("Position weight matrix:\n")
  print(round(object@pwm,4))
  cat("\n\nInformation content:\n")
  print(round(object@ic,4))
  cat("\n\nConsensus sequence:\n")
  print(object@consensus)
})

setMethod("print", signature(x="pwm"), function(x,...){
  print(round(x@pwm,4))
})

setMethod("show", "pwm", function(object){
  show(round(object@pwm,4))
})
