data_generation<-function(N, N_clusters, model, dataset1, dataset2){
    
    #arguments: -N: number of observations
    #           -N_clusters: number of Clusters
    #           -model: the data generating process
    #           -dataset1: randomization_data
    #           -dataset2: road_data
    #returns:    -simulation_data: dataframe 
    
    #merge data
    merged_data<-merge(dataset1, dataset2, by="desaid")
    merged_data<-merged_data[!is.na(merged_data$lndiffeall4mainancil),]
    merged_data$desaid<-as.numeric(merged_data$desaid)

    merged_data<-as.data.frame(merged_data)

    merged_data<-rename(merged_data, c("zdistancekec"="distance", "zkadesedyears"="head_education_years", "zpop"="population",
                         "zpercentpoorpra"="percentage_poor", "zkadesbengkoktotal"="head_salary", "zkadesage"="head_age",
                          "podeszhill"="mountainous", "totalmesjid"="mosques", "audit.x"="audit", 
                         "lndiffeall4mainancil"="ln_diff_items", "z4RABnumsubproj"="subprojects"))

    X_cont=as.data.frame(cbind(merged_data$population, merged_data$mosques, merged_data$totalallocation,
                               merged_data$subprojects, merged_data$percentage_poor,
                               merged_data$distance, merged_data$head_education_years, merged_data$head_age,
                               merged_data$head_salary))
    
    #drawing continuous variables

    mu<-colMeans(X_cont, na.rm=TRUE)

    cov<-var(X_cont, na.rm=TRUE)

    X_sim<-as.data.frame(rtmvnorm(n=N, mean = mu, sigma = cov, lower=rep(0,9), upper=c(rep(Inf,4), 1, rep(Inf,4))))    

    #logit
    
    glm.fit<-glm(mountainous ~ population+mosques+totalallocation+subprojects+percentage_poor+distance+head_education_years+head_age+head_salary, data = merged_data, family = "binomial")

    colnames(X_sim) <- c("population", "mosques", "totalallocation", "subprojects", "percentage_poor", "distance", "head_education_years", "head_age", "head_salary")


    mountainous<-c()
    X_sim$mountainous_prob<-predict(glm.fit,newdata=X_sim, type="response")
    for (i in 1:nrow(X_sim)){
        X_sim$mountainous[i]<-rbinom(n=1,size=1,prob=X_sim$mountainous_prob[i])
        }


    #clustering and randomization

    similarity_vars<-c("population", "mosques", "totalallocation", "percentage_poor", "mountainous")

    x_similarity<-as.data.frame(scale(X_sim[similarity_vars]))

    for (i in 1:nrow(x_similarity)){
        x_similarity$similarity_score[i]<-sum(x_similarity[i,])+rnorm(n=1, sd=1, mean=0)
    }

    X_sim<-cbind(X_sim, "similscore"=x_similarity$similarity_score)
    simulation_data<-X_sim[order(X_sim$similscore),]
    
    #specifying cluster size

    mean_cluster_size<-N/N_clusters

    min_cluster<-3/4*mean_cluster_size
    max_cluster<-5/4*mean_cluster_size

    cluster_vector<-rep(c(min_cluster, mean_cluster_size, max_cluster), times=N_clusters/3)

    simulation_data$subdistrict<-rep(1:N_clusters, times=c(cluster_vector))

    randomization_subdistricts<-rbinom(size=1, n=N_clusters, 0.48)

    cluster<-c()
    for (i in 1:N_clusters){
        cluster_temporary<-rep(randomization_subdistricts[i], times=cluster_vector[i])
        cluster<-c(cluster, cluster_temporary)   
    }
    simulation_data$audit<-cluster
    
    #specifying some further cluster variables

    idiosyncratic_main<-rnorm(mean=-0.05, sd=1, n=N_clusters)
    idiosyncratic_treat<-rnorm(mean=-0.05, sd=1, n=N_clusters)

    village_idiosyncratic_main<-c()

    for (i in 1:N_clusters){
        idiosyncratic_temporary<-rep(idiosyncratic_main[i], times=cluster_vector[i])
        village_idiosyncratic_main<-c(village_idiosyncratic_main, idiosyncratic_temporary)   
    }

    simulation_data$idiosyncratic_main<-village_idiosyncratic_main
    
      village_idiosyncratic_treat<-c()

    for (i in 1:N_clusters){
        idiosyncratic_temporary<-rep(idiosyncratic_treat[i], times=cluster_vector[i])
        village_idiosyncratic_treat<-c(village_idiosyncratic_treat, idiosyncratic_temporary)   
        }
    
    
    simulation_data$idiosyncratic_treat<-village_idiosyncratic_treat

    radio<-rtmvnorm(n=N_clusters, mean=1000, sigma=150000, lower=0)

    village_radio<-c()

    for (i in 1:N_clusters){
        radio_temporary<-rep(radio[i], times=cluster_vector[i])
        village_radio<-c(village_radio, radio_temporary)   
    }

    simulation_data$radio<-village_radio



    literacy<-rtmvnorm(n=N_clusters, mean=50000, sigma=300000000, lower=0)


    village_literacy<-c()

    for (i in 1:N_clusters){
        literacy_temporary<-rep(literacy[i], times=cluster_vector[i])
        village_literacy<-c(village_literacy, literacy_temporary)   
    }

    simulation_data$literacy<-village_literacy


    legalspending<-rtmvnorm(n=N_clusters, mean=50000, sigma=300000000, lower=0)

        village_legalspending<-c()

    for (i in 1:N_clusters){
        legalspending_temporary<-rep(legalspending[i], times=cluster_vector[i])
        village_legalspending<-c(village_legalspending, legalspending_temporary)   
    }

    simulation_data$legalspending<-village_legalspending


    #specying the covariance matrix structure

    covariance_matrix<-matrix(0, nrow=N, ncol=N)

        for (i in 1:N){
            for (j in ((i-max_cluster):(i+max_cluster))){
                if (j>0&j<=N){
                    if (simulation_data$subdistrict[i]==simulation_data$subdistrict[j]){
                    covariance_matrix[i,j]=0.2*0.001
                }
            }       
        }
        covariance_matrix[i,i]<-0.001
    }

    simulation_data$eps<-as.vector(rmvnorm(mean=rep(0, times=N), sigma=covariance_matrix, n=1))

    if (model=='simple')  {  

        simulation_data$ln_diff_items<-0.22-0.07*simulation_data$audit+
            0.0001*simulation_data$totalallocation-
            0.017*simulation_data$head_salary-
            0.016*simulation_data$audit*simulation_data$head_salary+
            simulation_data$eps

            simulation_data$tau<- -0.07-0.016*simulation_data$head_salary

        }

    if (model=='idiosyncratic') {



        simulation_data$ln_diff_items<-0.22-0.07*simulation_data$audit+
            0.0001*simulation_data$totalallocation-
            0.017*simulation_data$head_salary-
            0.016*simulation_data$audit*simulation_data$head_salary+
            simulation_data$idiosyncratic_main+
            simulation_data$idiosyncratic_treat*simulation_data$audit+
            simulation_data$eps

            simulation_data$tau<- -0.07-0.016*simulation_data$head_salary


            }

    if (model=='medium') {

        simulation_data$ln_diff_items<-0.27-0.11*simulation_data$audit+
            0.0001*simulation_data$totalallocation-
            0.0016*simulation_data$head_education_years+
            0.007*simulation_data$head_salary-
            0.0002*simulation_data$head_salary^2-
            0.01*simulation_data$audit*simulation_data$head_salary-
            0.001*simulation_data$audit*simulation_data$totalallocation+
            0.00003*simulation_data$head_salary*simulation_data$totalallocation+
            0.0002*simulation_data$audit*simulation_data$totalallocation*simulation_data$head_salary+
            0.0002*simulation_data$audit*simulation_data$head_salary^2+
            simulation_data$eps

            simulation_data$tau<- -0.11-0.01*simulation_data$head_salary-
            0.001*simulation_data$totalallocation+0.0002*simulation_data$totalallocation*
            simulation_data$head_salary+0.0002*simulation_data$head_salary^2


        }


    if (model=='complex') {

        simulation_data$ln_diff_items<-0.5-
        0.003*simulation_data$head_age-
        0.05*simulation_data$audit-
        0.01*simulation_data$percentage_poor-
        0.001*simulation_data$population-
        0.01*simulation_data$mosques-
        0.0005*simulation_data$totalallocation+
        0.005*simulation_data$head_salary-
        0.001*simulation_data$distance-
        0.002*simulation_data$subprojects+
        0.0002*simulation_data$audit*simulation_data$distance-
        0.2*simulation_data$audit*simulation_data$percentage_poor-
        0.02*simulation_data$audit*simulation_data$population-
        0.004*simulation_data$audit*simulation_data$mosques-
        0.002*simulation_data$audit*simulation_data$totalallocation-
        0.006*simulation_data$audit*simulation_data$head_salary+
        0.008*simulation_data$audit*simulation_data$subprojects+
        simulation_data$eps

        simulation_data$tau<- -0.05+
        0.0002*simulation_data$distance-
        0.2*simulation_data$percentage_poor-
        0.02*simulation_data$population-
        0.004*simulation_data$mosques-
        0.002*simulation_data$totalallocation-
        0.006*simulation_data$head_salary+
        0.008*simulation_data$subprojects

        }   

    return(simulation_data)
    }
    

MSE_function<-function(n_test, n_train, N_clusters, n.trees, model_true){
    
    #arguments: -n_test: number of test observations
    #           -_train: number of training observations
    #           -N_clusters: number of clusters
    #           -n.trees: number of trees
    #           -model_true: the dgp model
    #returns:   -list of MSEs, coverage probabilities and dataframe for plotting
    
    #preparing lists
    
    MSE_honest<-list()
    MSE_adaptive<-list()
    MSE_linear<-list()
    coverage_honest<-list()
    coverage_adaptive<-list()
    
    #drawing training and test data
    
    simulation_data<-data_generation(N=n_train, N_clusters=N_clusters, model=model_true, dataset1=randomization_data, dataset2=road_data)
    test_data<-data_generation(N=n_test, N_clusters=N_clusters, model=model_true, dataset1=randomization_data, dataset2=road_data)
    
    #specifying number of variables used if dgp is idiosyncratic
    if (model_true=='idiosyncratic'){
        X_train=cbind(simulation_data$totalallocation, simulation_data$head_salary, simulation_data$head_education_years, simulation_data$mosques,
                      simulation_data$population, simulation_data$subprojects, simulation_data$percentage_poor, simulation_data$distance, simulation_data$head_age, simulation_data$radio, simulation_data$literacy, simulation_data$legalspending)
        X_test=cbind(test_data$totalallocation, test_data$head_salary, test_data$head_education_years, test_data$mosques,
                      test_data$population, test_data$subprojects, test_data$percentage_poor, test_data$distance, test_data$head_age, test_data$radio, test_data$literacy, test_data$legalspending)
            }

    else{
        X_train=cbind(simulation_data$totalallocation, simulation_data$head_salary, simulation_data$head_education_years, simulation_data$mosques,
                  simulation_data$population, simulation_data$subprojects, simulation_data$percentage_poor, simulation_data$distance, simulation_data$head_age)
        X_test<-cbind(test_data$totalallocation, test_data$head_salary, test_data$head_education_years, test_data$mosques,
                  test_data$population, test_data$subprojects, test_data$percentage_poor, test_data$distance, 
                  test_data$head_age)
    }
    
    #looping over estimated models
    
    
    for (model_est in c('simple', 'medium_misspecified', 'complex', 'medium')){
    
        for (cluster in c('No', 'Yes')){

        if (cluster=='Yes')    { 
            train_forest_honest<-causal_forest(X_train, simulation_data$ln_diff_items, simulation_data$audit, W.hat=0.48, clusters=simulation_data$subdistrict, mtry=7, min.node.size=5, equalize.cluster.weights=TRUE, num.trees=n.trees, seed=1)
            train_forest_adaptive<-causal_forest(X_train, simulation_data$ln_diff_items, simulation_data$audit, W.hat=0.48, mtry=7, clusters=simulation_data$subdistrict, min.node.size=5, equalize.cluster.weights=TRUE, honesty=FALSE, num.trees=n.trees, seed=1)
            }

        if (cluster=='No') {
            train_forest_honest<-causal_forest(X_train, simulation_data$ln_diff_items, simulation_data$audit, W.hat=0.48,  mtry=7, min.node.size=5, num.trees=n.trees, seed=1)
            train_forest_adaptive<-causal_forest(X_train, simulation_data$ln_diff_items, simulation_data$audit, W.hat=0.48, mtry=7,  min.node.size=5, honesty=FALSE, num.trees=n.trees, seed=1)
            }



        if (model_est=='simple'){

            linear_model<-lm(ln_diff_items~totalallocation+audit*head_salary, 
                             data=simulation_data)

            coefficients<-linear_model$coefficients

            tauhat_linear<-coefficients[3]+coefficients[5]*test_data$head_salary

            }



        if (model_est=='medium'){

            linear_model<-lm(ln_diff_items~head_education_years+audit*totalallocation*head_salary+audit*I(head_salary^2), data=simulation_data)

            coefficients<-linear_model$coefficients

            tauhat_linear<-coefficients[3]+
            coefficients[7]*test_data$totalallocation+
            coefficients[8]*test_data$head_salary+
            coefficients[10]*test_data$head_salary^2+
            coefficients[11]*test_data$head_salary*test_data$totalallocation  
            }



            if (model_est=='complex'){

                linear_model<-lm(ln_diff_items~head_age+audit*distance+
                                 audit*percentage_poor+audit*population+
                                 audit*mosques+audit*totalallocation+
                                 audit*head_salary+audit*subprojects, data=simulation_data)

                coefficients<-linear_model$coefficients

                tauhat_linear<-coefficients[3]+
                                coefficients[11]*test_data$distance+
                    coefficients[12]*test_data$percentage_poor+
                    coefficients[13]*test_data$population+
                    coefficients[14]*test_data$mosques+
                    coefficients[15]*test_data$totalallocation+
                    coefficients[16]*test_data$head_salary+
                    coefficients[17]*test_data$subprojects

                    }






              if (model_est=='medium_misspecified'){

                linear_model<-lm(ln_diff_items~audit*head_salary+audit*totalallocation+audit*I(head_salary^2)+audit*
                             I(totalallocation^2),
                             data=simulation_data)

                coefficients<-linear_model$coefficients

                tauhat_linear<-coefficients[2]+
                                coefficients[7]*test_data$head_salary+
                                coefficients[8]*test_data$totalallocation+
                                coefficients[9]*(test_data$head_salary)^2+
                                coefficients[10]*(test_data$totalallocation)^2

                    }

    
    #predictions
    honest_predict<-predict(train_forest_honest, X_test, estimate.variance=TRUE)
    adaptive_predict<-predict(train_forest_adaptive, X_test, estimate.variance=TRUE)

    tauhat_honest<-honest_predict$predictions
    tauhat_adaptive<-adaptive_predict$predictions

    variance_honest<-honest_predict$variance.estimates
    variance_adaptive<-adaptive_predict$variance.estimates
            
    #Confidence intervals

    CI_honest_lower<-tauhat_honest-qnorm(0.95, mean=0, sd=1)*sqrt(variance_honest)
    CI_honest_upper<-tauhat_honest+qnorm(0.95, mean=0, sd=1)*sqrt(variance_honest)
    CI_adaptive_lower<-tauhat_adaptive-qnorm(0.95, mean=0, sd=1)*sqrt(variance_adaptive)
    CI_adaptive_upper<-tauhat_adaptive+qnorm(0.95, mean=0, sd=1)*sqrt(variance_adaptive)

    coverage_honest[[cluster]]<-mean(CI_honest_lower<=test_data$tau&test_data$tau<=CI_honest_upper)
    coverage_adaptive[[cluster]]<-mean(CI_adaptive_lower<=test_data$tau&test_data$tau<=CI_adaptive_upper)


    MSE_linear[[model_est]]<-mean((test_data$tau-tauhat_linear)^2)
    MSE_honest[[cluster]]<-mean((test_data$tau-tauhat_honest)^2)
    MSE_adaptive[[cluster]]<-mean((test_data$tau-tauhat_adaptive)^2)
        
    if (model_est=='medium_misspecified'){
            tauhat_linear_misspecified<-tauhat_linear
        }

    if (model_est=='medium'){
        complete_data<-cbind(test_data, tauhat_honest, tauhat_linear, tauhat_linear_misspecified)
        }
    }
    }


return(list("MSE Low Complexity Estimated Linear Model"=MSE_linear[['simple']],
            "MSE Misspecified Medium Complexity Estimated Linear Model"=MSE_linear[['medium_misspecified']],
            "MSE Medium Complexity Estimated Linear Model"=MSE_linear[['medium']],
            "MSE High Complexity Estimated Linear Model"=MSE_linear[['complex']],
            "MSE Honest, Cluster-Robust CF"=MSE_honest[['Yes']],
            "MSE Honest, not Cluster-Robust CF"=MSE_honest[['No']],
            "MSE Adaptive, Cluster-Robust CF"=MSE_adaptive[['Yes']],
            "Coverage Probability Honest, Cluster-Robust CF"=coverage_honest[['Yes']],
            "Coverage Probability Honest, not Cluster-Robust CF"=coverage_honest[['No']],
            complete_data=complete_data))
}

table_function<-function(observations, n_test, models, n.trees, N_clusters){

mylist <- list()
MSE_df <- data.frame()
this<-matrix()
for (n in observations){
for (model in models){
    mylist[[sprintf("True model: %s. %s", model, n)]]<-MSE_function(n_test=n_test, n_train=n, n.trees=n.trees, model_true=model, N_clusters=N_clusters)[1:9]
    MSE_df <- do.call("cbind",mylist)
}
}
    return(MSE_df)
}

plot_coverage_function<-function(n_test, n_train, n.trees, trees, observations){
dataf<-data.frame(
num_trees=NA,
cluster=NA,
nocluster=NA)
#coverage_prob<-list()
for (i in index(trees)){
    element<-trees[i]
    coverage<-MSE_function(n_test=n_test, n_train=n_train, n.trees=element, 
                           model_true='simple', N_clusters=150)[c("Coverage Probability Honest, Cluster-Robust CF", 
                                                  "Coverage Probability Honest, not Cluster-Robust CF")]
    #coverage_prob[[i]]<-rbind(element, "cluster"=coverage[1],"nocluster"=coverage[2])
    #dataf <- do.call("cbind", data.frame(coverage_prob))
    dataf[i,"num_trees"]<-element
    dataf[i,"cluster"]<-coverage[1]
    dataf[i,"no_cluster"]<-coverage[2]
    }

p1<-ggplot(dataf, aes(trees)) + 
  geom_line(aes(y = cluster, color = "cluster"))+
  geom_line(aes(y=no_cluster, color = "no_cluster"))+
  geom_point(aes(y = cluster, color = "cluster"))+
  geom_point(aes(y=no_cluster, color = "no_cluster"))+
  xlab('Log of number of trees')+
  ylab('Coverage Probability')+
  labs(x = 'Log of number of trees',
         y = 'Coverage Probability',
         color = "Legend")+
  scale_x_continuous(breaks = c(trees), trans = scales::log_trans())+
  geom_hline(yintercept= 0.9, size=0.3, linetype = "dashed", color='red')+
  scale_color_manual(values = c('cluster' = 'turquoise',
                                  'no_cluster' = 'black'))+
    ylim(0.6, 1)+
    ggtitle('Coverage: Different Number of Trees')+
    theme(plot.title = element_text(size = 10))
    
legend<-get_legend(p1)

p1<- p1 + theme(legend.position="none")


dataf<-data.frame(
num_observations=NA,
cluster=NA,
nocluster=NA)
#coverage_prob<-list()
for (i in index(observations)){
    element<-observations[i]
    coverage<-MSE_function(n_test=n_test, n_train=element, n.trees=n.trees, model_true='simple', N_clusters=150)[c("Coverage Probability Honest, Cluster-Robust CF", 
                                                  "Coverage Probability Honest, not Cluster-Robust CF")]
    #coverage_prob[[i]]<-rbind(element, "cluster"=coverage[1],"nocluster"=coverage[2])
    #dataf <- do.call("cbind", data.frame(coverage_prob))
    dataf[i,"num_observations"]<-element
    dataf[i,"cluster"]<-coverage[1]
    dataf[i,"no_cluster"]<-coverage[2]
    }

p2<-ggplot(dataf, aes(observations)) + 
  geom_line(aes(y = cluster, colour = "cluster"))+
  geom_line(aes(y = no_cluster, colour= "no_cluster"))+
  geom_point(aes(y = cluster, colour = "cluster"))+
  geom_point(aes(y = no_cluster, colour= "no_cluster"))+
  xlab('Log of number of observations')+
  ylab('Coverage Probability')+
  scale_x_continuous(breaks = c(observations), trans = scales::log_trans())+
  theme(legend.position="none")+
  scale_colour_manual(values = c('cluster' = 'turquoise',
                                  'no_cluster' = 'black'))+
  geom_hline(yintercept=0.9, size=0.3, linetype="dashed", 
                color = "red")+
    ylim(0.6, 1)+
    ggtitle('Coverage: Different Number of Observations')+
    theme(plot.title = element_text(size = 10))

grid.arrange(p1, p2, legend, ncol=2, nrow = 2, 
             layout_matrix = rbind(c(1,2), c(3,3)),
             widths = c(2.7, 2.7), heights = c(2.5, 0.5))

}

figure2<-function(n_test, n_train, n.trees=n.trees, model, N_clusters){

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
    
cf_empirical_function<-function(dataset1, dataset2){

merged_data<-merge(dataset1, dataset2, by="desaid")
merged_data<-merged_data[!is.na(merged_data$lndiffeall4mainancil),]
merged_data$desaid<-as.numeric(merged_data$desaid)

merged_data<-as.data.frame(merged_data)

merged_data<-rename(merged_data, c("zdistancekec"="distance", "zkadesedyears"="head_education_years", "zpop"="population",
                     "zpercentpoorpra"="percentage_poor", "zkadesbengkoktotal"="head_salary", "zkadesage"="head_age",
                      "podeszhill"="mountainous", "totalmesjid"="mosques", "audit.x"="audit", "kecnum.x"="subdistrict", 
                     "lndiffeall4mainancil"="ln_diff_items", "z4RABnumsubproj"="subprojects"))
covariates<-c('distance', 'head_education_years', 'population', 'percentage_poor', 'head_salary',
              'head_age', 'mountainous', 'mosques', 'subprojects')

X=merged_data[,covariates]
cf_honest_clustered<-causal_forest(X=X, Y=merged_data$ln_diff_items, W=merged_data$audit, clusters=merged_data$subdistricts, tune.parameters=mtry, seed=1)
cf_honest_not_clustered<-causal_forest(X=X, Y=merged_data$ln_diff_items, W=merged_data$audit, tune.parameters=mtry, seed=1)
cf_adaptive<-causal_forest(X=X, Y=merged_data$ln_diff_items, W=merged_data$audit, honesty=FALSE, tune.parameteters=mtry, seed=1)
    
variable_importance_honest_clustered=variable_importance(cf_honest_clustered)
variable_importance_adaptive=variable_importance(cf_adaptive)
#tau_hat_clustered = predict(cf_clustered)$predictions
#hist(tau_hat_clustered)
#tau_hat_not_clustered = predict(cf_not_clustered)$predictions
#hist(tau_hat_not_clustered)

#sd(tau_hat_clustered)
#sd(tau_hat_not_clustered)

ATE_honest_clustered<-average_treatment_effect(cf_honest_clustered)
ATE_honest_not_clustered<-average_treatment_effect(cf_honest_not_clustered)
ATE_adaptive<-average_treatment_effect(cf_adaptive)

tauhat_honest_clustered<-predict(cf_honest_clustered)$predictions
tauhat_honest_not_clustered<-predict(cf_honest_not_clustered)$predictions
tauhat_adaptive<-predict(cf_adaptive)$predictions

return(list(ATE_honest_clustered=ATE_honest_clustered, ATE_honest_not_clustered=ATE_honest_not_clustered,
            ATE_adaptive=ATE_adaptive, tauhat_honest_clustered=tauhat_honest_clustered, 
           tauhat_honest_not_clustered=tauhat_honest_not_clustered, 
           tauhat_adaptive=tauhat_adaptive, cf_honest_clustered=cf_honest_clustered, 
           cf_honest_not_clustered=cf_honest_not_clustered, cf_adaptive=cf_adaptive, X=X, 
           merged_data=merged_data, 
           variable_importance_honest_clustered=variable_importance_honest_clustered,
           variable_importance_adaptive=variable_importance_adaptive))
}

ATE_function<-function(dataset1, dataset2){
empirical_results<-cf_empirical_function(randomization_data, road_data)
df<-data.frame()

df['ATE','Honest cluster-robust CF']<-empirical_results$ATE_honest_clustered[1]
df['Standard Error','Honest cluster-robust CF']<-empirical_results$ATE_honest_clustered[2]
df['ATE','Honest non-cluster-robust CF']<-empirical_results$ATE_honest_not_clustered[1]
df['Standard Error','Honest non-cluster-robust CF']<-empirical_results$ATE_honest_not_clustered[2]
df['ATE', 'Adaptive non-cluster-robust CF']<-empirical_results$ATE_adaptive[1]
df['Standard Error', 'Adaptive non-cluster-robust CF']<-empirical_results$ATE_adaptive[2]
    
    return(df)
}

boxplots<-function(forest, variable_importance){
empirical_results=cf_empirical_function(randomization_data, road_data)
cf=empirical_results[forest]
X=empirical_results$X
merged_data=empirical_results$merged_data
variable_importance=unlist(empirical_results[variable_importance], use.names=FALSE)
important_variables=which(variable_importance > mean(variable_importance))

if (forest=='cf_honest_clustered'){

cf_imp<-causal_forest(X=X[,important_variables], Y=merged_data$ln_diff_items,
                      W=merged_data$audit, clusters=merged_data$subdistricts)
    }

if (forest=='cf_adaptive'){

cf_imp<-causal_forest(X=X[,important_variables], Y=merged_data$ln_diff_items,
                      W=merged_data$audit, tune.parameters="all")
    }

tau_hat_clustered<-predict(cf_imp)$predictions

plot_data<-as.data.frame(cbind(X[,important_variables], tau_hat_clustered))
names<-names(plot_data[,-(ncol(plot_data))])
names
p<-list()
for (i in index(names)){
    name<-names[i]
    p[[i]]<-ggplot(plot_data, aes_string(x=plot_data[,name], y=tau_hat_clustered)) + 
          geom_boxplot(aes(group = cut_width(plot_data[,name], ((min(plot_data[,name]))+(max(plot_data[,name])))/10)))+
            geom_smooth(se=FALSE)+
    xlab(name)+
    ylab('CATE')
    }
do.call(grid.arrange, c(p, ncol=2))
    }