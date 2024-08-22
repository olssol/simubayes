server <- function(input, output, session) {

    source("des_tools.R")

    y = reactive({
        if(input$borrow == 'control arm' & input$outcome =='binary'){
            return(as.numeric(unlist(strsplit(input$y,","))))
        }})
    N = reactive({
        if(input$borrow == 'control arm' & input$outcome =='binary'){
            return(as.numeric(unlist(strsplit(input$N,","))))
        }})

    y.m = reactive({
        if(input$borrow == 'control arm' & input$outcome =='normal'){
            return(as.numeric(unlist(strsplit(input$y.m,","))))
        }})

    y.se = reactive({
        if(input$borrow == 'control arm' & input$outcome =='normal'){
            return(as.numeric(unlist(strsplit(input$y.se,","))))
        }})

    y.dif = reactive({
        if(input$borrow == 'treatment difference between arms'){
            return(as.numeric(unlist(strsplit(input$y.dif,","))))
        }})
    y.dif.se = reactive({
        if(input$borrow == 'treatment difference between arms'){
            return(as.numeric(unlist(strsplit(input$y.dif.se,","))))
        }})
    logOR = reactive({
        if(input$borrow == 'log OR'){
            return(as.numeric(unlist(strsplit(input$logOR,","))))
        }})
    logOR.se = reactive({
        if(input$borrow == 'log OR'){
            return(as.numeric(unlist(strsplit(input$logOR.se,","))))
        }})

    rt1 = reactive({as.numeric(unlist(strsplit(input$rt1,",")))})
    rt2_temp = reactive({as.numeric(unlist(strsplit(input$rt2,",")))})
    rt2 = reactive({c(rt2_temp()[1], diff(rt2_temp()))})

    success = reactive({as.numeric(unlist(strsplit(input$success,",")))})
    futility = reactive({as.numeric(unlist(strsplit(input$futility,",")))})

    prior.hist = eventReactive(input$start, {
        if(input$prior=='Meta-Analytic Predictive Prior'){
            if(input$borrow=='treatment difference between arms'){
                priors = stb_tl_bayes_para_rmap(dtype = input$outcome,
                                                borrow = 'diff',
                                                y.dif = y.dif(),
                                                y.dif.se = y.dif.se(),
                                                v = input$v,
                                                m = input$m,
                                                s = input$s,
                                                weight_noninfo = input$w,
                                                prior1_noninfo = input$p1.noninfo,
                                                prior2_noninfo = input$p2.noninfo,
                                                cores = 1)
            }
            if(input$borrow=='control arm'){
                if(input$outcome =='normal'){
                    priors = stb_tl_bayes_para_rmap(dtype = input$outcome,
                                                    borrow = 'control',
                                                    y.m = y.m(),
                                                    y.se = y.se(),
                                                    v = input$v,
                                                    m = input$m,
                                                    s = input$s,
                                                    weight_noninfo = input$w,
                                                    prior1_noninfo = input$p1.noninfo,
                                                    prior2_noninfo = input$p2.noninfo,
                                                    cores = 1)
                }
                if(input$outcome =='binary'){
                    priors = stb_tl_bayes_para_rmap(dtype = input$outcome,
                                                    borrow = 'control',
                                                    y = y(),
                                                    N = N(),
                                                    v = input$v,
                                                    m = input$m,
                                                    s = input$s,
                                                    weight_noninfo = input$w,
                                                    prior1_noninfo = input$p1.noninfo,
                                                    prior2_noninfo = input$p2.noninfo,
                                                    cores = 1)
                }
            }
            if(input$borrow=='log OR'){
                priors = stb_tl_bayes_para_rmap(dtype = input$outcome,
                                                borrow = 'logOR',
                                                logOR = logOR(),
                                                logOR.se = logOR.se(),
                                                v = input$v,
                                                m = input$m,
                                                s = input$s,
                                                weight_noninfo = input$w,
                                                prior1_noninfo = input$p1.noninfo,
                                                prior2_noninfo = input$p2.noninfo,
                                                cores = 1)
            }
        }
        if(input$prior=='Power Prior'){
            if(input$borrow=='treatment difference between arms'){
                priors = stb_tl_bayes_para_pp(y.dif = y.dif(),
                                              y.dif.se = y.dif.se(),
                                              dtype = input$outcome,
                                              borrow = 'diff',
                                              prior1 = input$a,
                                              prior2 = input$b,
                                              power.param = input$power,
                                              cores=1)
            }
            if(input$borrow=='control arm'){
                if(input$outcome == 'normal'){
                    priors = stb_tl_bayes_para_pp(y.m = y.m(),
                                                  y.se = y.se(),
                                                  dtype = input$outcome,
                                                  borrow = 'control',
                                                  prior1 = input$a,
                                                  prior2 = input$b,
                                                  power.param = input$power,
                                                  cores=1)
                }
                if(input$outcome == 'binary'){
                    priors = stb_tl_bayes_para_pp(y = y(),
                                                  N = N(),
                                                  dtype = input$outcome,
                                                  borrow = 'control',
                                                  prior1 = input$a,
                                                  prior2 = input$b,
                                                  power.param = input$power,
                                                  cores=1)
                }
            }
            if(input$borrow=='log OR'){
                priors = stb_tl_bayes_para_pp(logOR = logOR(),
                                              logOR.se = logOR.se(),
                                              dtype = input$outcome,
                                              borrow = 'logOR',
                                              prior1 = input$a,
                                              prior2 = input$b,
                                              power.param = input$power,
                                              cores=1)
            }
        }
        return(priors)
    })

    lst_para = eventReactive(prior.hist(), {
        if(input$outcome == 'normal'){
            if(input$borrow == 'treatment difference between arms'){
                lst = list(sample_size     = 30,
                           ratio_by_arm    = rt1(),
                           ratio_by_stage  = rt2(),
                           prior_by_arm    = list(prior.hist()),
                           m_s_control     = c(input$mu_c, input$sd_c),
                           m_s_treat       = c(input$mu_t, input$sd_t),
                           n_post          = 10000,
                           test_direction  = input$direct,
                           decision_h0     = input$target,
                           decision_suc    = success(),
                           decision_fut    = futility(),
                           borrow          = 'diff')
            }
            if(input$borrow == 'control arm'){
                lst = list(sample_size     = 30,
                           ratio_by_arm    = rt1(),
                           ratio_by_stage  = rt2(),
                           prior_by_arm    = list(prior.hist(), c(input$t.a, input$t.b)),
                           m_s_control     = c(input$mu_c, input$sd_c),
                           m_s_treat       = c(input$mu_t, input$sd_t),
                           n_post          = 10000,
                           test_direction  = input$direct,
                           decision_h0     = input$target,
                           decision_suc    = success(),
                           decision_fut    = futility(),
                           borrow          = 'control')
            }
        }

        if(input$outcome == 'binary'){
            if(input$borrow == 'treatment difference between arms'){
                lst = list(sample_size     = 30,
                           ratio_by_arm    = rt1(),
                           ratio_by_stage  = rt2(),
                           p_by_arm        = c(input$pi.c, input$pi.t),
                           prior_by_arm    = list(prior.hist()),
                           n_post          = 10000,
                           test_direction  = input$direct,
                           decision_h0     = input$target,
                           decision_suc    = success(),
                           decision_fut    = futility(),
                           borrow          = 'diff',
                           test_statistic  = 'diff')
            }
            if(input$borrow == 'log OR'){
                lst = list(sample_size     = 30,
                           ratio_by_arm    = rt1(),
                           ratio_by_stage  = rt2(),
                           p_by_arm        = c(input$pi.c, input$pi.t),
                           prior_by_arm    = list(prior.hist()),
                           n_post          = 10000,
                           test_direction  = input$direct,
                           decision_h0     = input$target,
                           decision_suc    = success(),
                           decision_fut    = futility(),
                           borrow          = 'logOR',
                           test_statistic  = 'logOR')
            }
            if(input$borrow == 'control arm'){
                lst = list(sample_size     = 30,
                           ratio_by_arm    = rt1(),
                           ratio_by_stage  = rt2(),
                           p_by_arm        = c(input$pi.c, input$pi.t),
                           prior_by_arm    = list(prior.hist(), c(input$t.a, input$t.b)),
                           n_post          = 10000,
                           test_direction  = input$direct,
                           decision_h0     = input$target,
                           decision_suc    = success(),
                           decision_fut    = futility(),
                           borrow          = 'control',
                           test_statistic  = input$teststat)
            }
        }
        return(lst)
    })
                                        #
    ss.seq = reactive({seq(input$ssc1, input$ssc2, input$SSi)})

    octab = eventReactive(lst_para(),{
        stb_tl_get_oc_shiny(lst_design = lst_para(), dtype = input$outcome, ss_seq = ss.seq(), n_rep = input$ntrials)
    })

    lst_para_h0 = eventReactive(prior.hist(), {
        if(input$outcome == 'normal'){
            if(input$borrow == 'treatment difference between arms'){
                lst = list(sample_size     = 30,
                           ratio_by_arm    = rt1(),
                           ratio_by_stage  = rt2(),
                           prior_by_arm    = list(prior.hist()),
                           m_s_control     = c(input$mu_c, input$sd_c),
                           m_s_treat       = c(input$mu_c+ input$target, input$sd_t),
                           n_post          = 10000,
                           test_direction  = input$direct,
                           decision_h0     = input$target,
                           decision_suc    = success(),
                           decision_fut    = futility(),
                           borrow          = 'diff')
            }
            if(input$borrow == 'control arm'){
                lst = list(sample_size     = 30,
                           ratio_by_arm    = rt1(),
                           ratio_by_stage  = rt2(),
                           prior_by_arm    = list(prior.hist(), c(input$t.a, input$t.b)),
                           m_s_control     = c(input$mu_c, input$sd_c),
                           m_s_treat       = c(input$mu_c + input$target, input$sd_t),
                           n_post          = 10000,
                           test_direction  = input$direct,
                           decision_h0     = input$target,
                           decision_suc    = success(),
                           decision_fut    = futility(),
                           borrow          = 'control')
            }
        }

        if(input$outcome == 'binary'){
            if(input$borrow == 'treatment difference between arms'){
                lst = list(sample_size     = 30,
                           ratio_by_arm    = rt1(),
                           ratio_by_stage  = rt2(),
                           p_by_arm        = c(input$pi.c, input$pi.c + input$target),
                           prior_by_arm    = list(prior.hist()),
                           n_post          = 10000,
                           test_direction  = input$direct,
                           decision_h0     = input$target,
                           decision_suc    = success(),
                           decision_fut    = futility(),
                           borrow          = 'diff',
                           test_statistic  = 'diff')
            }
            if(input$borrow == 'log OR'){
                lst = list(sample_size     = 30,
                           ratio_by_arm    = rt1(),
                           ratio_by_stage  = rt2(),
                           p_by_arm        = c(input$pi.c, input$pi.c + input$target),
                           prior_by_arm    = list(prior.hist()),
                           n_post          = 10000,
                           test_direction  = input$direct,
                           decision_h0     = input$target,
                           decision_suc    = success(),
                           decision_fut    = futility(),
                           borrow          = 'logOR',
                           test_statistic  = 'logOR')
            }
            if(input$borrow == 'control arm'){
                lst = list(sample_size     = 30,
                           ratio_by_arm    = rt1(),
                           ratio_by_stage  = rt2(),
                           p_by_arm        = c(input$pi.c, input$pi.c + input$target),
                           prior_by_arm    = list(prior.hist(), c(input$t.a, input$t.b)),
                           n_post          = 10000,
                           test_direction  = input$direct,
                           decision_h0     = input$target,
                           decision_suc    = success(),
                           decision_fut    = futility(),
                           borrow          = 'control',
                           test_statistic  = input$teststat)
            }
        }
        return(lst)
    })


    output$OC_tab = renderTable(octab(), rownames = T, colnames = F)

    octab_h0 = eventReactive(lst_para_h0(),{
        stb_tl_get_oc_shiny(lst_design = lst_para_h0(), dtype = input$outcome, ss_seq = ss.seq(), n_rep = input$ntrials)
    })

    output$OC_plot = renderPlot({stb_tl_plot_oc_shiny(lst_para(), octab())})

    output$OC_plot_h0 = renderPlot({stb_tl_plot_type1_shiny(lst_para_h0(), octab_h0())})

    output$type1 = renderTable({data.frame(total_sample_size = octab_h0()[1,],
                                           type_1_error = octab_h0()[sapply(length(rt2()), function(i){paste('success rate stage', i)}),])})
    output$bayes1 = renderPlot({stb_tl_get_bayes_update_shiny(lst_para(), dtype = input$outcome, ss_seq = ss.seq())})

    session$onSessionEnded(function() {
                stopApp()
            })
}
