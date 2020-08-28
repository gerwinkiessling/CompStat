figure1<-function(n_test, n_train, n.trees=n.trees, model, N_clusters){

plot_data<-MSE_function(n_test=n_test, n_train=n_train, n.trees=n.trees, model_true=model, N_clusters=N_clusters)$complete_data
plot_data<-plot_data[!(plot_data$head_salary>14| plot_data$totalallocation>180),]


minimum<-min(plot_data$tau, plot_data$tauhat_honest, plot_data$tauhat_linear, plot_data$tauhat_linear_misspecified)
maximum<-max(plot_data$tau, plot_data$tauhat_honest, plot_data$tauhat_linear, plot_data$tauhat_linear_misspecified)

p1<-ggplot(plot_data, aes(head_salary, totalallocation))+
    geom_point(aes(colour = tau)) +
    scale_colour_gradientn(colours = terrain.colors(10),limits=c(minimum, maximum))+
    ggtitle(TeX('Medium Complexity Model: True Effect $\\tau(x)$'))+
    theme(plot.title = element_text(size = 10))+
    theme(legend.title=element_blank())



legend <- get_legend(p1)

p1<- p1 + theme(legend.position="none")

p2<-ggplot(plot_data, aes(head_salary, totalallocation))+
    geom_point(aes(colour = tauhat_honest)) +
    scale_colour_gradientn(colours = terrain.colors(10), limits=c(minimum,maximum))+
    theme(legend.position="none")+
    ggtitle('Estimate from Honest Cluster-Robust\n Causal Forest')+
    theme(plot.title = element_text(size = 10))

p3<-ggplot(plot_data, aes(head_salary, totalallocation))+
    geom_point(aes(colour = tauhat_linear)) +
    scale_colour_gradientn(colours = terrain.colors(10), limits=c(minimum,maximum))+
    theme(legend.position="none")+
    ggtitle('Estimate from Medium Complexity \n Linear Model')+
    theme(plot.title = element_text(size = 10))

p4<-ggplot(plot_data, aes(head_salary, totalallocation))+
    geom_point(aes(colour = tauhat_linear_misspecified)) +
    scale_colour_gradientn(colours = terrain.colors(10), limits=c(minimum,maximum))+
    theme(legend.position="none")+
    ggtitle('Estimate from Misspecified Medium \n Complexity Linear Model')+
    theme(plot.title = element_text(size = 10))


grid.arrange(p1, p2, p3, p4, legend, ncol=3, nrow = 2, 
             layout_matrix = cbind(c(1,2), c(3,4), c(5,5)),
             widths = c(2.5, 2.5, 0.5), heights = c(2.5,2.5))
    }