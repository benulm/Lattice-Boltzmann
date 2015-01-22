a <- read.csv("~/ETH/Master/7. Semester/Lattice Boltzi/exercises/Lattice-Boltzmann/output/1u_profile.dat", header=FALSE)
b <- read.csv("~/ETH/Master/7. Semester/Lattice Boltzi/exercises/Lattice-Boltzmann/output/8001u_profile.dat", header=FALSE)
c <- read.csv("~/ETH/Master/7. Semester/Lattice Boltzi/exercises/Lattice-Boltzmann/output/20001u_profile.dat", header=FALSE)
d <- read.csv("~/ETH/Master/7. Semester/Lattice Boltzi/exercises/Lattice-Boltzmann/output/30001u_profile.dat", header=FALSE)
e <- read.csv("~/ETH/Master/7. Semester/Lattice Boltzi/exercises/Lattice-Boltzmann/output/50001u_profile.dat", header=FALSE)
f <- read.csv("~/ETH/Master/7. Semester/Lattice Boltzi/exercises/Lattice-Boltzmann/output/80001u_profile.dat", header=FALSE)
g <- read.csv("~/ETH/Master/7. Semester/Lattice Boltzi/exercises/Lattice-Boltzmann/output/100001u_profile.dat", header=FALSE)
h <- read.csv("~/ETH/Master/7. Semester/Lattice Boltzi/exercises/Lattice-Boltzmann/output/130001u_profile.dat", header=FALSE)
i <- read.csv("~/ETH/Master/7. Semester/Lattice Boltzi/exercises/Lattice-Boltzmann/output/200001u_profile.dat", header=FALSE)
j <- read.csv("~/ETH/Master/7. Semester/Lattice Boltzi/exercises/Lattice-Boltzmann/output/400001u_profile.dat", header=FALSE)


#matplot((a[,1])^(1/7),(a[,1]), type="l", col="green")
matplot(a[,2],a[,1], type="l", col="red", xlim=c(0,0.1))
lines(b[,2],b[,1], type="l", col="green")
lines(c[,2],c[,1], type="l", col="cyan")
lines(d[,2],d[,1], type="l",col="blue")
lines(e[,2],e[,1], type="l")
lines(f[,2],f[,1], type="l", col="red")
lines(h[,2],h[,1], type="l", col="orange")
lines(i[,2],i[,1], type="l", col="cyan")
lines(j[,2],j[,1], type="l", col="green")


u_h <- read.csv("~/ETH/Master/7. Semester/Lattice Boltzi/exercises/Lattice-Boltzmann/output/u_h.dat", header=FALSE)
matplot(u_h[,1],u_h[,2], type="l", col="red")


a <- read.csv("~/ETH/Master/7. Semester/Lattice Boltzi/exercises/Lattice-Boltzmann/output/1I_profile.dat", header=FALSE)
b <- read.csv("~/ETH/Master/7. Semester/Lattice Boltzi/exercises/Lattice-Boltzmann/output/8001I_profile.dat", header=FALSE)
c <- read.csv("~/ETH/Master/7. Semester/Lattice Boltzi/exercises/Lattice-Boltzmann/output/20001I_profile.dat", header=FALSE)
d <- read.csv("~/ETH/Master/7. Semester/Lattice Boltzi/exercises/Lattice-Boltzmann/output/30001I_profile.dat", header=FALSE)
e <- read.csv("~/ETH/Master/7. Semester/Lattice Boltzi/exercises/Lattice-Boltzmann/output/50001I_profile.dat", header=FALSE)
f <- read.csv("~/ETH/Master/7. Semester/Lattice Boltzi/exercises/Lattice-Boltzmann/output/80001I_profile.dat", header=FALSE)
g <- read.csv("~/ETH/Master/7. Semester/Lattice Boltzi/exercises/Lattice-Boltzmann/output/100001I_profile.dat", header=FALSE)
h <- read.csv("~/ETH/Master/7. Semester/Lattice Boltzi/exercises/Lattice-Boltzmann/output/130001I_profile.dat", header=FALSE)
i <- read.csv("~/ETH/Master/7. Semester/Lattice Boltzi/exercises/Lattice-Boltzmann/output/200001I_profile.dat", header=FALSE)
j <- read.csv("~/ETH/Master/7. Semester/Lattice Boltzi/exercises/Lattice-Boltzmann/output/400001I_profile.dat", header=FALSE)

#matplot((a[,1])^(1/7),(a[,1]), type="l", col="green")
matplot(a[,2],a[,1], type="l", col="red", xlim=c(0,1))
lines(b[,2],b[,1], type="l", col="green")
lines(c[,2],c[,1], type="l", col="cyan")
lines(d[,2],d[,1], type="l",col="blue")
lines(e[,2],e[,1], type="l")
lines(f[,2],f[,1], type="l", col="red")
lines(g[,2],g[,1], type="l", col="purple")
lines(h[,2],h[,1], type="l", col="orange")
lines(i[,2],i[,1], type="l", col="cyan")
lines(j[,2],j[,1], type="l", col="green")








