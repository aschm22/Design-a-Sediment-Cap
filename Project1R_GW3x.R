rm(list=ls())

Cs = 129500000                     # sediment solid concentration, ug/kg (given)
logKoc = 5.2                       # organic-water partition coefficient (given)
foc = 0.05                         # organic content of the sediment (given)
Cpw = Cs/(foc*10^logKoc)           # porewater concentration in the sediment
remedialsed = 7500                 # remedial goal of sediment ug/kg
remedialpw = remedialsed / (foc* 10^logKoc)

Dm = 6.54e-6 * 3600 * 24           # this converts molecular dispersivity into cm2/day
Dbiopw = 100 / 365                 # this converts Dbiopw into cm2/day
Dbios = 1 / 365                    # converts Dbios into cm2/day
vgw = 3*198 / 365                    # converts groundwater velocity into cm/day
alpha = 0.1                        # hydraulic dispersivity (given)
theta_bio = 0.4                    # bioturbation layer porosity (given)
rho_bio = 1.56                     # bioturbation layer material density, kg/L (given)
theta_cil = 0.5                    # chemical isolation layer porosity 
rho_cil = 1.25                     # chemical isolation layer material density, kg/L

dx = 1                             # model grid length, cm
nbio = 15                          # total grids for bioturbation layer sim
ncil = 15                          # total grids for chemical isolation layer sim
tt = 365 * 250                     # total simulation time (250 years) in days

vbio = vgw / theta_bio             # velocity of porewater in bioturbation layer
vcil = vgw / theta_cil             # velocity of porewater in chemical isolation layer
Dtbiopw = Dm + Dbiopw + alpha * vbio  # total porewater effective diffusion coefficient in bioturbation layer
Dtbios = Dbios                     # total solid effective diffusion coefficient in bioturbation layer
Dtcilpw = Dm + alpha * vcil        # total porewater effective diffusion coefficient in chemical isolation layer
logKd = 1.9                        # adsorption coefficient
Kd = 10^logKd                      # log transformation
logKf = 7.21                       # activated carbon Freundlich
Kf = 10^logKf                      # log transformation
N = 0.82                           # Freundlich 1/n

dt1bio = dx / vbio                 # time step within bioturbation layer
dt2bio = dx * dx / Dtbiopw / 2
dtbio = min(dt1bio, dt2bio)

dt1cil = dx / vcil                 # time step within chemical isolation layer
dt2cil = dx * dx / Dtcilpw / 2
dtcil = min(dt1cil, dt2cil)

dt = min(dtbio, dtcil) / 2
nt = ceiling(tt/dt)
dt = tt / nt

C_cil = array(1e-15, ncil+1)       # initial conditions
C_bio = array(0, nbio)
C_bio_s = C_bio
C_bio_old = C_bio
C_output <- C_cil

fac = 0.0060                       # activated carbon fraction (play around with this)

plot(-30:-1, append(C_cil[2:(ncil+1)], C_bio), ylim=c(1e-15, 1e3), log='y', xlab='Depth (cm)', ylab='Cw (ug/L)', type = 'l', col='blue')
abline(v = -15, lty = 3, col='gray', lwd = 2)           # interface between bioturbation and chemical isolation layer
abline(h = remedialpw, lty = 3, col='red', lwd = 2)
t.prt <- 0
col.plot <- 0

while (C_cil[16] < remedialpw) {
  
  
  
  for (k in 1:nt) {
    t.prt = t.prt + dt
    C_cil[1] <- Cpw
    for (i in 2:ncil) {
      R = 1 + Kd * rho_cil * (1-fac) / theta_cil + Kf * N * C_cil[i]^(N-1) * rho_cil * fac / theta_cil
      A <- vcil * dt / R / dx
      B <- Dtcilpw * dt / R / dx / dx
      C_cil[i] <- (1-A-2*B) * C_cil[i] + (A + B) * C_cil[i-1] + B*C_cil[i+1]
    }
    R = 1 + Kd * rho_bio * (1-fac) / theta_cil + Kf * N * C_cil[ncil + 1]^(N-1) * rho_cil * fac / theta_cil
    A <- vcil * dt / R / dx
    B <- Dtcilpw * dt / R / dx / dx
    C_cil[ncil+1] <- (1 - A - B) * C_cil[ncil+1] + (A+B) * C_cil[ncil]
    
    R = 1 + Kd * rho_bio / theta_bio
    A <- vbio * dt / R / dx
    B <- Dtbiopw * dt / R / dx / dx
    C_bio_old = C_bio
    
    C_bio[1] <- (1-A-B) * C_bio[1] + A * C_cil[ncil+1] + B * C_bio[2]
    for (i in 2:(nbio-1)) {
      C_bio[i] <- (1-A-2*B) * C_bio[i] + (A+B) * C_bio[i-1] + B * C_bio[i+1]
    }
    C_bio[nbio] <- (1-A-B) * C_bio[nbio] + (A+B) * C_bio[nbio-1]
    
    Ss <- rho_bio * Kd * (C_bio - C_bio_old)
    B <- Dtbios * dt / dx / dx
    C_bio_s[1] <- (1-B) * C_bio_s[1] + B * C_bio_s[2] + Ss[1]
    for (i in 2:(nbio-1)) {
      C_bio_s[i] <- (1-2*B) * C_bio_s[i] + B * C_bio_s[i-1] + C_bio_s[i+1] + Ss[i]
    }
    C_bio_s[nbio] <- (1-B) * C_bio_s[nbio] + B * C_bio_s[nbio-1] + Ss[nbio]
    C_bio = C_bio_s / Kd
    
    if (t.prt >= 3650) {
      lines(-30:-1, append(C_cil[2:(ncil+1)], C_bio), col=col.plot)
      t.prt = 0
      col.plot <- col.plot+1
    }
  }
  
  fac = fac - 0.0001
}
