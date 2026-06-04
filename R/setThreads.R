#' @title Set the number of OpenMP threads
#' @description Controls the number of threads used by the C++ backend. Simple 
#' package internal to access RcppParallel function, `setThreadOptions`.
#' @param numThreads Integer. Number of threads to use for task scheduling. 
#' Call defaultNumThreads() to determine the the default value used for "auto".
#' @param stackSize Integer. Stack size (in bytes) to use for worker threads. 
#' The default used for "auto" is 2MB on 32-bit systems and 4MB on 64-bit 
#' systems (note that this parameter has no effect on Windows).
#' @importFrom RcppParallel setThreadOptions
#' @export
setOpenMPThreads <- function(numThreads = "auto", stackSize = "auto") {
  RcppParallel::setThreadOptions(numThreads = numThreads, stackSize = stackSize)
}