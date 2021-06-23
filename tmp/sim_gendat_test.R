set.seed(124)

data_I <- gendat_HeCE(l3n=5,l2n=5,l1n=7,
                      mu0=3,mu1=0.01,
                      sg00=0.05,sg11=0.005,sg01=-0.005, 
                      sp0=0.2,
                      sp=0.2,
                      se=0.2,
                      delta1=-0.08)

ggplot(data_I) + 
  geom_point(aes(time,p, color=factor(person)))+
  facet_grid(rows=vars(group)) +
  geom_line(aes(time,p, color=factor(person)))+
  geom_point(aes(time,g))

ggplot(data_I) + 
  geom_point(aes(time,y, color=factor(person)))+
  facet_grid(rows=vars(group)) +
  geom_line(aes(time,y, color=factor(person)))+
  geom_point(aes(time,g))

data_IGP <- gendat_HeCEGP(l3n=5,l2n=5,l1n=3,
                      mu0=3,mu1=0.01,
                      #sg00=0.05,sg11=0.005,sg01=-0.005, 
                      sp0=0.01,
                      sp=0.2,
                      se=0.01,
                      sigma=0.2,
                      corr=0.05,
                      delta1=-0.8,
                      delta2 = 0.8)

ggplot(data_IGP) + 
  geom_point(aes(time,p, color=factor(person)))+
  facet_grid(rows=vars(group)) +
  geom_line(aes(time,p, color=factor(person)))+
  geom_point(aes(time,g))

ggplot(data_IGP) + 
  geom_point(aes(time,y, color=factor(person)))+
  facet_grid(rows=vars(group)) +
  geom_line(aes(time,y, color=factor(person)))+
  geom_point(aes(time,g))

data_II <- gendat_HoCE(l3n=5,l2n=5,l1n=3,
                      mu0=3,mu1=0.01,
                      sg00=0.05,sg11=0.005,sg01=-0.005, 
                      sp0=0.01,
                      sp=0.2,
                      se=0.01,
                      delta1=-0.8)

ggplot(data_II) + 
  geom_point(aes(time,p, color=factor(person)))+
  facet_grid(rows=vars(group)) +
  geom_line(aes(time,p, color=factor(person)))+
  geom_point(aes(time,g))

ggplot(data_II) + 
  geom_point(aes(time,y, color=factor(person)))+
  facet_grid(rows=vars(group)) +
  geom_line(aes(time,y, color=factor(person)))+
  geom_point(aes(time,g))

data_IIGP <- gendat_HoCEGP(l3n=5,l2n=5,l1n=3,
                           mu0=3,mu1=0.01,
                           #sg00=0.05,sg11=0.005,sg01=-0.005, 
                           sp0=0.01,
                           sp=0.2,
                           se=0.01,
                           sigma=0.2,
                           corr=0.05,
                           delta1=-0.8,
                           delta2 = 0.8)

ggplot(data_IIGP) + 
  geom_point(aes(time,p, color=factor(person)))+
  facet_grid(rows=vars(group)) +
  geom_line(aes(time,p, color=factor(person)))+
  geom_point(aes(time,g))

ggplot(data_IIGP) + 
  geom_point(aes(time,y, color=factor(person)))+
  facet_grid(rows=vars(group)) +
  geom_line(aes(time,y, color=factor(person)))+
  geom_point(aes(time,g))
