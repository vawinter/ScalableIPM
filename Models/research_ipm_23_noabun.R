ipm <- nimbleCode({
  for (i in 1:hwb.N) {
    HWB[i] ~ dbern(hwb.p[i])
    logit(hwb.p[i]) <- hwb.beta1 + hwb.beta2 * hwb.Year2020[i] + 
      hwb.beta3 * hwb.Year2021[i] + hwb.beta4 * hwb.Year2022[i] + 
      hwb.beta5 * hwb.Year2023[i] + hwb.beta6 * hwb.doy.scale[i] + 
      hwb.beta7 * hwb.doy.2[i] + hwb.u[hwb.wmu[i]]
  }
  hwb.beta1 ~ dnorm(0, sd = 1)
  hwb.beta2 ~ dnorm(0, sd = 1)
  hwb.beta3 ~ dnorm(0, sd = 1)
  hwb.beta4 ~ dnorm(0, sd = 1)
  hwb.beta5 ~ dnorm(0, sd = 1)
  hwb.beta6 ~ dnorm(0, sd = 1)
  hwb.beta7 ~ dnorm(0, sd = 1)
  hwb.sigma ~ dunif(0, 10)
  for (j in 1:hwb.J) {
    hwb.u[j] ~ dnorm(0, sd = hwb.sigma)
  }
  for (t in 1:Nyears) {
    for (u in 1:female.n.wmu) {
      aug31.hwb[t, u] <- expit(hwb.beta1 + hwb.beta2 * 
                                 (t == 1) + hwb.beta3 * (t == 2) + hwb.beta4 * 
                                 (t == 3) + hwb.beta5 * (t == 4) + hwb.beta6 * 
                                 hwb.aug31 + hwb.beta7 * hwb.aug31.2 + hwb.u[u])
    }
  }
  for (k in 1:ph.N) {
    PHratio[k] ~ dgamma(shape = ph.alpha[k], scale = ph.theta[k])
    ph.alpha[k] <- ph.mu[k] * ph.disp
    ph.theta[k] <- 1/ph.disp
    log(ph.mu[k]) <- ph.beta1 + ph.beta2 * ph.Year2020[k] + 
      ph.beta3 * ph.Year2021[k] + ph.beta4 * ph.Year2022[k] + 
      ph.beta5 * ph.Year2023[k] + ph.beta6 * ph.doy.scale[k] + 
      ph.beta7 * ph.doy.2[k] + ph.u[ph.wmu[k]]
  }
  ph.beta1 ~ dnorm(0, sd = 1)
  ph.beta2 ~ dnorm(0, sd = 1)
  ph.beta3 ~ dnorm(0, sd = 1)
  ph.beta4 ~ dnorm(0, sd = 1)
  ph.beta5 ~ dnorm(0, sd = 1)
  ph.beta6 ~ dnorm(0, sd = 1)
  ph.beta7 ~ dnorm(0, sd = 1)
  ph.disp ~ dunif(0, 10)
  ph.sigma.u ~ dunif(0, 1)
  for (j in 1:ph.J) {
    ph.u[j] ~ dnorm(0, sd = ph.sigma.u)
  }
  for (t in 1:Nyears) {
    for (u in 1:female.n.wmu) {
      aug31.ppb[t, u] <- exp(ph.beta1 + ph.beta2 * (t == 
                                                      1) + ph.beta3 * (t == 2) + ph.beta4 * (t == 3) + 
                               ph.beta5 * (t == 4) + ph.beta6 * ppb.aug31 + 
                               ph.beta7 * ppb.aug31.2 + ph.u[u])
    }
  }
  telem.beta.int[1] ~ dnorm(0, sd = 0.5)
  telem.beta.age[1] ~ dnorm(0, sd = 0.5)
  telem.sigma ~ dunif(0, 5)
  for (w in 1:female.telem.wmu) {
    telem.beta.wmu[w] ~ dnorm(0, sd = telem.sigma)
  }
  telem.month.sigma ~ dunif(0, 5)
  for (m in 1:12) {
    telem.beta.month[m] ~ dnorm(0, sd = telem.month.sigma)
  }
  for (i in 1:telem.nind) {
    for (y in telem.year.start[i]:telem.year.end[i]) {
      for (m in telem.first[i]:telem.last[i]) {
        s.kf[i, y, m] <- telem.beta.int[1] + telem.beta.age[1] * 
          telem.juvenile[i, y, m] + inprod(telem.beta.wmu[1:female.telem.wmu], 
                                           telem.wmu[i, 1:female.telem.wmu]) + telem.beta.month[m]
        status[i, y, m] ~ dbern(prob = icloglog(s.kf[i, 
                                                     y, m]))
      }
    }
  }
  for (m in 1:12) {
    storage[1, 1, m] <- icloglog(telem.beta.int[1] + telem.beta.wmu[1] + 
                                   telem.beta.month[m])
    storage[1, 2, m] <- icloglog(telem.beta.int[1] + telem.beta.age[1] + 
                                   telem.beta.wmu[1] + telem.beta.month[m])
  }
  for (j in 2:(female.telem.wmu)) {
    for (m in 1:12) {
      storage[j, 1, m] <- icloglog(telem.beta.int[1] + 
                                     telem.beta.wmu[j] + telem.beta.month[m])
      storage[j, 2, m] <- icloglog(telem.beta.int[1] + 
                                     telem.beta.age[1] + telem.beta.wmu[j] + telem.beta.month[m])
    }
  }
  for (j in 1:female.telem.wmu) {
    avg.ad.s.kf[j] <- prod(storage[j, 1, 1:12])
  }
  for (j in 1:female.telem.wmu) {
    juvenile_part1[j] <- storage[j, 2, 1] * storage[j, 2, 
                                                    1]
    juvenile_part2[j] <- prod(storage[j, 2, 1:5])
    adult_part[j] <- prod(storage[j, 1, 6:10])
    avg.juv.s.kf[j] <- (juvenile_part1[j] * juvenile_part2[j]) * 
      adult_part[j]
  }
  for (j in 1:male.n.wmu) {
    juvenile_male_1[j] <- storage[j, 2, 1] * storage[j, 2, 
                                                     1]
    juvenile_male_2[j] <- prod(storage[j, 2, 1:4])
    juv.male.adj[j] <- (juvenile_male_1[j] * juvenile_male_2[j])
  }
  for (t in 1:(female.n.occasions - 1)) {
    for (u in 1:female.n.wmu) {
      female.h.ad.wmu[t, u] ~ dbeta(shape1 = 2, shape2 = 50)
      female.h.juv.wmu[t, u] ~ dbeta(shape1 = 2, shape2 = 50)
    }
  }
  male.rrate.j ~ dnorm(0.71, sd = 0.072)
  male.rrate.a ~ dnorm(0.87, sd = 0.039)
  male.seber.recov[1] ~ dunif(0, 1)
  male.seber.recov[2] ~ dunif(0, 1)
  male.juvenile.effect ~ dnorm(0, sd = 0.5)
  for (t in 1:(male.n.occasions - 1)) {
    male.time.effect[t] ~ dnorm(0, sd = 0.5)
  }
  male.sigma ~ dunif(0, 10)
  for (i in 1:10) {
    male.wmu.effect[i] ~ dnorm(0, sd = male.sigma)
  }
  for (i in 1:male.nind) {
    for (t in male.f[i]:(male.n.occasions - 1)) {
      logit(male.s[i, t]) <- male.juvenile.effect * male.I[i, 
                                                           t] + inprod(male.time.effect[1:4], male.time.param[t, 
                                                                                                              1:4]) + male.wmu.effect[male.wmu[i]]
      male.r[i, t] <- male.seber.recov[1] * male.I[i, t] * 
        male.rrate.j * male.II[i] + male.seber.recov[1] * 
        male.I[i, t] * (1 - male.II[i]) + male.seber.recov[2] * 
        (1 - male.I[i, t]) * male.rrate.a * male.II[i] + 
        male.seber.recov[2] * (1 - male.I[i, t]) * (1 - 
                                                      male.II[i])
    }
  }
  for (t in 1:(male.n.occasions - 1)) {
    for (u in 1:male.n.wmu) {
      logit(male.s.juv.wmu[t, u]) <- male.juvenile.effect + 
        inprod(male.time.effect[1:4], male.time.param[t, 
                                                      1:4]) + male.wmu.effect[u]
      logit(male.s.ad.wmu[t, u]) <- inprod(male.time.effect[1:4], 
                                           male.time.param[t, 1:4]) + male.wmu.effect[u]
    }
  }
  for (t in 1:(male.n.occasions - 1)) {
    for (u in 1:male.n.wmu) {
      male.h.juv.wmu[t, u] <- (1 - male.s.juv.wmu[t, u]) * 
        male.seber.recov[1]
      male.h.ad.wmu[t, u] <- (1 - male.s.ad.wmu[t, u]) * 
        male.seber.recov[2]
    }
  }
  for (i in 1:male.nind) {
    male.z[i, male.f[i]] <- 1
    for (t in (male.f[i] + 1):male.n.occasions) {
      male.z[i, t] ~ dbern(male.mu1[i, t])
      male.mu1[i, t] <- male.s[i, t - 1] * male.z[i, t - 
                                                    1]
      male.y[i, t] ~ dbern(male.mu2[i, t])
      male.mu2[i, t] <- male.r[i, t - 1] * (male.z[i, t - 
                                                     1] - male.z[i, t])
    }
  }
  # for (u in 1:male.n.wmu) {
  #   N.lambda.ad.male[u] <- (th.year1.male.ad[u])/male.h.ad.wmu[1, 
  #                                                              u]
  #   male.N.ad[1, u] ~ dpois(N.lambda.ad.male[u])
  #   N.lambda.juv.male[u] <- (th.year1.male.juv[u])/male.h.juv.wmu[1, 
  #                                                                 u]
  #   male.N.juv[1, u] ~ dpois(N.lambda.juv.male[u])
  #   male.recruitment[1, u] <- recruitment[1, u] * juv.male.adj[u]
  # }
  # for (t in 2:Nyears) {
  #   for (u in 1:male.n.wmu) {
  #     male.N.ad.Survived[t, u] ~ dbin(prob = male.s.ad.wmu[t - 
  #                                                            1, u], size = male.N.ad[t - 1, u])
  #     male.N.juv.Survived[t, u] ~ dbin(prob = male.s.juv.wmu[t - 
  #                                                              1, u], size = male.N.juv[t - 1, u])
  #     male.N.ad[t, u] <- (male.N.ad.Survived[t, u] + male.N.juv.Survived[t, 
  #                                                                        u])
  #     male.recruitment[t, u] <- recruitment[t, u] * juv.male.adj[u]
  #     male.N.juv[t, u] ~ dpois(male.recruitment[t - 1, 
  #                                               u])
  #     harvest.ad.spring[t, u] ~ dbin(prob = male.h.ad.wmu[t, 
  #                                                         u], size = male.N.ad[t, u])
  #     harvest.juv.spring[t, u] ~ dbin(prob = male.h.juv.wmu[t, 
  #                                                           u], size = male.N.juv[t, u])
  #   }
  # }
  # for (u in 1:female.n.wmu) {
  #   N.lambda.ad.female[u] <- (th.year1.female.ad[u])/female.h.ad.wmu[1, 
  #                                                                    u]
  #   female.N.ad[1, u] ~ dpois(N.lambda.ad.female[u])
  #   N.lambda.juv.female[u] <- (th.year1.female.juv[u])/female.h.juv.wmu[1, 
  #                                                                       u]
  #   female.N.juv[1, u] ~ dpois(N.lambda.juv.female[u])
  #   recruitment[1, u] <- ((female.N.ad[1, u] * aug31.hwb[1, 
  #                                                        u]) * aug31.ppb[1, u])/2
  # }
  # for (t in 2:Nyears) {
  #   for (u in 1:female.n.wmu) {
  #     recruitment[t, u] <- ((female.N.ad[t, u] * aug31.hwb[t, 
  #                                                          u]) * aug31.ppb[t, u])/2
  #     female.N.ad.Survived[t, u] ~ dbin(size = female.N.ad[t - 
  #                                                            1, u], prob = avg.ad.s.kf[u])
  #     female.N.juv.Survived[t, u] ~ dbin(size = female.N.juv[t - 
  #                                                              1, u], prob = avg.juv.s.kf[u])
  #     female.N.ad[t, u] <- (female.N.ad.Survived[t, u] + 
  #                             female.N.juv.Survived[t, u])
  #     female.N.juv[t, u] ~ dpois(recruitment[t, u])
  #     harvest.ad.fall[t, u] ~ dbin(prob = female.h.ad.wmu[t, 
  #                                                         u], size = female.N.ad[t, u])
  #     harvest.juv.fall[t, u] ~ dbin(prob = female.h.juv.wmu[t, 
  #                                                           u], size = female.N.juv[t, u])
  #   }
  # }
})
