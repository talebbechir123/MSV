library(deSolve)

# ODE function for 4 variants, with competition form b_{i,j} = r_j / r_i
# and mutation m_{i,j} depending on distances in the VPg region
virus_ode <- function(t, state, params) {
  # state = (V1, V2, V3, V4)
  V1 <- state[1]
  V2 <- state[2]
  V3 <- state[3]
  V4 <- state[4]
  
  # Unpack parameters
  r1  <- params["r1"]; r2  <- params["r2"]
  r3  <- params["r3"]; r4  <- params["r4"]
  K   <- params["K"]
  mut <- params["mut"]  # single-mutation rate
  
  # For convenience, define each V to avoid rewriting
  V <- c(V1, V2, V3, V4)
  r <- c(r1, r2, r3, r4)
  
  # We assume 2-locus differences:
  # Distances: (1) NN->DN or NH, (2) DN->DH or NN, etc...
  # We'll define a function dist(i,j) that returns 0 if i=j, 1 if differ by 1 mut, else 2
  dist <- matrix(c(
    0,1,1,2,
    1,0,2,1,
    1,2,0,1,
    2,1,1,0
  ), nrow=4, byrow=TRUE)
  
  # Define mutation rates m_{i,j}
  # If distance = 1 => m
  # If distance = 2 => m^2
  # If distance = 0 => 0
  m_ij <- matrix(0, nrow=4, ncol=4)
  for(i in 1:4) {
    for(j in 1:4) {
      if(dist[i,j] == 1) {
        m_ij[i,j] <- mut
      } else if(dist[i,j] == 2) {
        m_ij[i,j] <- mut^2
      } else {
        m_ij[i,j] <- 0
      }
    }
  }
  
  # We also define b_{i,j} = r_j / r_i
  # We'll store them in a matrix b[i,j]
  b <- matrix(0, nrow=4, ncol=4)
  for(i in 1:4) {
    for(j in 1:4) {
      if(i != j) {
        b[i,j] <- r[j] / r[i]
      }
    }
  }
  
  # Now compute dV_i/dt
  dV <- numeric(4)
  sumV <- sum(V)
  
  for(i in 1:4) {
    # competition term: r_i * V_i(t) * ( 1 - [V_i + sum_{j != i} b_{i,j} V_j ] / K )
    comp_part <- V[i] * r[i] * (1 - (V[i] + sum(b[i, -i]*V[-i]))/K )
    
    # mutation term: sum_{j!=i} m_{i,j} * [ V_j - V_i ]
    # expand that to sum_{j} because if j=i, it's zero anyway
    mut_part <- 0
    for(j in 1:4) {
      mut_part <- mut_part + m_ij[i,j]*(V[j] - V[i])
    }
    
    dV[i] <- comp_part + mut_part
  }
  
  list(dV)
}

#### Set up parameters
params <- c(
  r1=1.05,  # Example: variant 1's growth rate
  r2=0.91,  # variant 2
  r3=1.04,  # variant 3
  r4=0.99,  # variant 4
  K=1e6,    # carrying capacity
  mut=1e-5  # single point mutation probability
)

#### Initial conditions (based on the inoculum frequencies)
# Suppose the inoculum frequencies were NN=0.32, DN=0.22, NH=0.22, DH=0.24
# We'll scale by total virus = 100 for the initial time
V_tot0 <- 100
y0 <- c(V_tot0*0.32, V_tot0*0.22, V_tot0*0.22, V_tot0*0.24)

#### Solve over 0..35 days
times <- seq(0, 35, by=0.5)

out <- ode(y=y0, times=times, func=virus_ode, parms=params)

# Convert to data frame
out <- as.data.frame(out)
colnames(out) <- c("time","V1","V2","V3","V4")

# Compute frequencies
freq <- out
freq[,2:5] <- freq[,2:5] / rowSums(freq[,2:5])

# Plot frequencies over time
matplot(freq$time, freq[,2:5], type="l", lty=1, lwd=2,
        xlab="Time (days)", ylab="Frequency",
        col=c("blue","red","green","purple"))
legend("topright", legend=c("NN","DN","NH","DH"),
       col=c("blue","red","green","purple"), lty=1, lwd=2, bty="n")


