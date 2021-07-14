## r(mat_t) plots for simulation 

true_rtmax <- rep(rep(c(1,0.92,0.88,0.83),each = 7),4)
ce_type <- rep(c("HetCE","HomCE","HetCE GP","HomCE GP"),each=28)
d1 <- rep(rep(c(0,-0.03,-0.05,-0.08),each=7),4)
Model <- rep(rep(c("HetCEM", "HomCEM", "GP", 
               "HetCEM GP", "HomCEM GP", "GP GP","CEM"),4),4)

rtmax <- c(
          # HetCE
          1,1,1,0.99,1,1,1,
          0.93,0.97,0.90,0.90,0.98,0.90,0.95,
          0.89,0.93,0.87,0.87,0.91,0.85,0.92,
          0.85,0.88,0.83,0.82,0.85,0.81,0.89,
          
          # HomCE
          1,1,1,0.99,1,1,1,
          1,0.94,1,0.98,0.91,1,1,
          0.99,0.89,0.99,0.97,0.86,1,0.99,
          0.98,0.84,0.98,0.96,0.80,1,0.99,
          
          # HetCE GP
          rep(NA,28),
          
          # HomCE GP
          1.01,1,1.03,1,1,1,1.03,
          1.01,0.96,1.02,0.99,0.94,0.99,1.02,
          1,0.92,1.02,0.99,0.9,0.98,1.02,
          1,0.86,1.01,0.98,0.85,0.96,1.02)

dat <- data.frame(ce_type,d1,true_rtmax,Model,rtmax)

ggplot(dat,aes(d1,rtmax))+
  facet_wrap(vars(ce_type))+
  geom_point(aes(color=Model)) +
  geom_line(aes(color=Model,linetype=Model)) +
  theme_bw() +
  geom_point(aes(d1,true_rtmax))+
  xlab(expression(delta [1])) +
  ylab("r(tmax)")

ggplot(dat,aes(d1,rtmax))+
  facet_wrap(vars(Model))+
  #facet_grid(rows = vars(group), cols = vars(person),labeller = label_both) +
  geom_point(aes(color=ce_type)) +
  geom_line(aes(color=ce_type)) +
  theme_bw() +
  geom_point(aes(d1,true_rtmax)) +
  xlab(expression(delta [1])) +
  ylab("r(tmax)")
