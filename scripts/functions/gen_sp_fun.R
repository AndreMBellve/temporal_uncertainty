#Notes on the function arguments: 

#"species.type" needs to be set to multiplicative to achieve distributions similar to those produced by MaxEnt. 

#"rescale" set to true rescales the probability of occurrence to be between 0 and 1 across the entire map. This might not be appropriate though as presumably in some time periods the habitat might not be favourable. Maybe it means it sets it at a maximum of 1 and a minimum of zero though?

#"rescale.each.response" forces the probability of occurrence to be between 0 and 1 for each time period, essentially forcing there to be optimal habitats between the environmental range observed. Setting this to FALSE essentially forces the niche to be constant through time.

#Altered function to work with terra for greater efficiency
gen_sp_fun <- function(env.rast, 
                       parameters,
                       rescale = TRUE, 
                       formula = NULL, 
                       species.type = "multiplicative", 
                       rescale.each.response = TRUE,
                       plot = FALSE) 
{
  message("Generating virtual species environmental suitability...\n")
  approach <- "response"
  

## Defensive programming ---------------------------------------------------

  #Checking object type is correct
  if (!(is(env.rast, "SpatRaster"))) {
    stop("env.rast must be a SpatRaster object")
  }
  
  #Setting min/max if they are incorrect - to be done in pre-processing of layers
  if (any(is.infinite(as.vector(minmax(env.rast))))) {
    "Print"
    env.rast <- setMinMax(env.rast)
  }
  
  #Checking the # of layers in SpatRaster object matches the number of functions...
  n.l <- nlyr(env.rast)
  if (n.l != length(parameters)) {
    stop("Provide as many layers in env.rast as functions on parameters")
  }
  
  #...and that the names match too
  if (any(!(names(parameters) %in% names(env.rast)) | 
          !(names(env.rast) %in% names(parameters)))) {
    stop("Layer names and names of parameters must be identical")
  }
  for (i in 1:length(parameters)) {
    if (any(!(c("fun", "args") %in% names(parameters[[i]])))) {
      stop("The structure of parameters does not seem correct. \n          Please provide function and arguments for variable '", 
           names(parameters)[i], "'. See help(generateSpFromFun) for more details.", 
           sep = "")
    }
    
    #Defensive checks
    test <- tryCatch(match.fun(parameters[[i]]$fun), error = function(c) "error")
    if (!inherits(test, "function")) {
      stop(paste("The function ", parameters[[i]]$fun, 
                 " does not exist, please verify spelling.", 
                 sep = ""))
    }
    
    
    if (any(!(names(parameters[[i]]$args) %in% names(formals(fun = test))))) {
      stop(paste("Arguments of variable '", names(parameters)[i], 
                 "' (", paste(names(parameters[[i]]$args), collapse = ", "), 
                 ") do not match arguments of the associated function\n\n                 List of possible arguments for this function: ", 
                 paste(names(formals(fun = test)), collapse = ", "), 
                 sep = ""))
    }
    rm(test)
  }
  
  #Rescaling the response variables and environmental suitability
  if (rescale.each.response) {
    message(" - The response to each variable was rescaled between 0 and 1. To\n            disable, set argument rescale.each.response = FALSE\n")
  }
  
  if (rescale) {
    message(" - The final environmental suitability was rescaled between 0 and 1.\n            To disable, set argument rescale = FALSE\n")
  }

  tic("suitability raster")
  
  #Creating suitability raster based on the desired functions
  suitab.raster <- sapply(names(env.rast), 
                          FUN = function(y) {
    app(env.rast[[y]], 
        
        fun = function(x) {
           
           do.call(match.fun(parameters[[y]]$fun), 
                   args = c(list(x),
                            parameters[[y]]$args))
           })
    }) %>% 
    rast()
  toc()
  
  for(var in names(env.rast)) {
    parameters[[var]]$min <- env.rast[[var]]@ptr$range_min
    parameters[[var]]$max <- env.rast[[var]]@ptr$range_max
  }
  
  
  #Updated, but untested
  if(rescale.each.response) {
    suitab.raster <- sapply(names(suitab.raster), 
                                  function(y){
                                    (suitab.raster[[y]] - suitab.raster[[y]]@ptr$range_min)/(suitab.raster[[y]]@ptr$range_max - 
                                                                                          suitab.raster[[y]]@ptr$range_min)
                                  }) %>% 
      rast()
  }
  
  #Combining raster results by either taking the product of the values or summing them
  if (is.null(formula)) {
    if (species.type == "multiplicative") {
      formula <- paste(names(suitab.raster), collapse = " * ")
      suitab.raster <- app(suitab.raster,
                            fun = prod)
    }
    else if (species.type == "additive") {
      formula <- paste(names(suitab.raster), collapse = " + ")
      suitab.raster <- app(suitab.raster,
                            fun = terra::sum)
    }
    else stop("If you do not provide a formula, please choose either species.type = 'additive' or 'multiplicative'")
  } else {
    if (any(!(all.vars(reformulate(formula)) %in% names(suitab.raster)))) {
      stop("Please verify that the variable names in your formula are correctly spelled")
    }
    else if (any(!(names(suitab.raster) %in% all.vars(reformulate(formula))))) {
      stop("Please verify that your formula contains all the variables of your input raster stack")
    }
    else {
      custom.fun <- NULL
      eval(parse(text = paste("custom.fun <- function(", 
                              paste(names(suitab.raster), collapse = ", "), 
                              ") {", formula, "}")))
      suitab.raster <- app(suitab.raster, 
                                       fun = custom.fun)
      print(formula)
    }
  }
  
  #Rescaling the values to be between 0 and 1 probabilities.
  if (rescale) {
    suitab.raster <- (suitab.raster - suitab.raster@ptr$range_min)/(suitab.raster@ptr$range_max - suitab.raster@ptr$range_min)
  }
  
  
  #Combining results to return list object
  results <- list(approach = approach, 
                  details = list(variables = names(parameters), 
                                 formula = formula, 
                                 rescale.each.response = rescale.each.response,
                                 rescale = rescale, 
                                 parameters = parameters), 
                  suitab.raster = suitab.raster)
  
  #Plotting?
  if (plot) {
    plot(results$suitab.raster, 
         main = "Environmental suitability of the virtual species")
  }
  
  #Returning results
  class(results) <- append("virtualspecies", class(results))
  return(results)
}

#Bench marking shows that the terra version is roughly 2x faster (564 v.s.287). Likely greater gains with more layers.
