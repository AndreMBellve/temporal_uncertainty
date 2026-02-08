#Script that reads in a list of the virtual species response functions identical to those created for our preservation biases study. However, this list is a subset of the species used for that study as we felt some species were redundant and these represented the environmental niche space well

vs_resp_ls <- list(
  #Defining imitation response curves (Blarina carolinensis)
  "bla_car" = formatFunctions(elevation = c(fun = 'dsigmoid',
                                            curve = -0.005,
                                            halfway = 200),
                              
                              precip_meanAnnual = c(fun = 'logisticFun',
                                                    alpha = -50,
                                                    beta = 1000),
                              
                              precip_seasonality = c(fun = 'dsigmoid',
                                                     curve = -20,
                                                     halfway = 125),
                              
                              slope = c(fun = 'dsigmoid',
                                        curve = -20,
                                        halfway = 45),
                              
                              temp_meanAnnual = c(fun = 'dnorm',
                                                  mean = 15,
                                                  sd = 7.5),
                              
                              temp_seasonality = c(fun = "logisticFun",
                                                   alpha = -30,
                                                   beta = 250)),
  
  #Defining imitation response curves (Lepus californicus)
  "lep_cal" = formatFunctions(elevation = c(fun = 'dsigmoid',
                                            curve = 20,
                                            halfway = -999),
                              
                              precip_meanAnnual = c(fun = 'dsigmoid',
                                                    curve = -0.004,
                                                    halfway = 1000),
                              
                              precip_seasonality = c(fun = 'logisticFun',
                                                     alpha = -2,
                                                     beta = 12.5),
                              
                              slope = c(fun = 'dsigmoid',
                                        curve = -20,
                                        halfway= 45),
                              
                              temp_meanAnnual = c(fun = 'dnorm',
                                                  mean = 17,
                                                  sd = 8),
                              
                              temp_seasonality = c(fun = "dgamma",
                                                   shape = 3.5,
                                                   rate = 0.005)),
  
  
  #Defining imitation response curves (Microtus pennsylvanicus)
  "mic_pen" = formatFunctions(elevation = c(fun = 'dgevd',
                                            location = 0,
                                            scale = 3000,
                                            shape = -0.4),
                              
                              precip_meanAnnual = c(fun = 'dweibull',
                                                    shape = 1.2,
                                                    scale = 2000),
                              
                              precip_seasonality = c(fun = 'inv_d_sigmoid',
                                                     curve = 25,
                                                     max_val = 165),
                              
                              slope = c(fun = 'dsigmoid',
                                        curve = -20,
                                        halfway = 45),
                              
                              temp_meanAnnual = c(fun = 'dgevd',
                                                  location = 5,
                                                  scale = 8,
                                                  shape = 0.55),
                              
                              temp_seasonality = c(fun = "inv_i_sigmoid",
                                                   curve = 100,
                                                   max_val = 2010)),
  
  
  #Defining imitation response curves (Neotoma albigula)
  "neo_alb" = formatFunctions(elevation = c(fun = 'logisticFun',
                                            alpha = -1000,
                                            beta = 400),
                              
                              precip_meanAnnual = c(fun = 'inv_d_sigmoid',
                                                    max_val = 5000,
                                                    curve = 80),
                              
                              precip_seasonality = c(fun = 'logisticFun',
                                                     alpha = -10,
                                                     beta = 30),
                              
                              slope = c(fun = 'dsigmoid',
                                        curve = -20,
                                        halfway= 45),
                              
                              temp_meanAnnual = c(fun = 'logisticFun',
                                                  alpha = -1.8,
                                                  beta = 17),
                              
                              temp_seasonality = c(fun = "logisticFun",
                                                   alpha = -100,
                                                   beta = 500)))