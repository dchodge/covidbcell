#' Calculate the trajectories for the validations data given the stanfit object for each individual
#' @param stanfit The stanfit object
#' @param foldername The foldername
cal_traj_val <- function(stanfit, foldername) {
    
    data_meta <- readRDS(here::here("outputs", "data_clean", foldername, "df_meta.RDS"))
    data_listA <- readRDS(here::here("outputs", "data_clean", foldername, "list_stan_data.RDS"))

    PIDs <- data_meta$PID 

    param_uncert <- stanfit$draws(
        c("theta_t",   
            "prop_asc", "deltaAg", "deltaP", "deltaL", "deltaA", "deltaG")) %>%
        spread_draws(theta_t[j, i, s], prop_asc, deltaAg, deltaP, deltaL, deltaA, deltaG) %>% 
        mutate(deltaM = 0.001, .before = "deltaP") %>%
        mutate(deltaAg = 1 / (deltaAg * 30), deltaP = 1 / deltaP, deltaL = 1 / deltaL, deltaA = 1 / deltaA, deltaG = 1 / deltaG) %>% 
        pivot_wider(names_from = "s", values_from = "theta_t") %>% 
        select(j:.draw, `1`:`4`, prop_asc, deltaAg:deltaG) %>% 
        rename(t = j, v = i) %>% rename(param1 = `1`, param2 = `2`, param3 = `3`, param4 = `4`)


    param_uncert_sum <- param_uncert %>% ungroup %>% summarise(across(param1:deltaG, ~mean(.x)), .by = c("t", "v")) 

    M <- 100
    df_uncert_all <- map_df(1:length(PIDs), 
        function(s) {
            init_i <- data_listA$y0_init[s, ]
            names(init_i) <- c("B", "P", "A1", "G", "L", "A2")

            param_uncert_sum_i <- param_uncert %>% filter(t == data_listA$x_time[s], v == data_listA$x_vac[s])

            df_uncert <- map_df(1:M,
                ~ode_results_vac(param_uncert_sum_i[.x, 6:16], 365, init_i) %>% as.data.frame  )  %>%
                select(time, B, P, A1, A2) %>% mutate(A = A1 + A2) %>% 
                mutate(i = map(1:M, ~rep(.x, 366)) %>% unlist) %>% 
                mutate(PID = PIDs[s]) %>% pivot_longer(B:A, names_to = "state", values_to = "value") %>%
                group_by(time, PID, state) %>% mutate(n = M)
        }
    )

    df_uncert_all_add <- df_uncert_all %>% left_join(data_meta)

}

#' Calculate the trajectories for the validations data given the stanfit object for each individual
run_val_trajectories <- function(stanfit, foldername) {
    if (foldername == "sid_wu_s") {
        antigen_str <- "Ancestral spike"
        antigen_file <- "s"
    } else {
        antigen_str <- "Ancestral RBD"
        antigen_file <- "rbd"
    }
    df_uncert_all_labs_s <- cal_traj_val(stanfit, foldername)
    saveRDS(df_uncert_all_labs_s, here::here("outputs", "stanfit",  paste0("traj_val_",  antigen_file, ".RDS")))
}

#' Calculate the poluation-level trajectories for the validations data given the stanfit object for each individual
#' 
plot_val_traj_pop <- function() {

     # Get the trajectories for each individual
    df_uncert_ind_s <- readRDS(here::here("outputs", "stanfit",  paste0("traj_val_",  "s", ".RDS")))  %>% mean_qi(value) %>% mutate(antigen = "Spike model")
    df_uncert_ind_rbd <- readRDS(here::here("outputs", "stanfit",  paste0("traj_val_",  "rbd", ".RDS")))  %>% mean_qi(value) %>% mutate(antigen = "RBD model")
    df_uncert_ind <- bind_rows(df_uncert_ind_s, df_uncert_ind_rbd) %>% filter(state %in% c("B", "P", "A")) %>%
        mutate(state_variable = recode(state, B = "Memory B-cell conc.", P = "Plasmablast conc.", A = "sVNT to ancestral variant")) %>% mutate(n = 4000)

    # Get the calibration data
    data_i_all_s <- get_data_ind("sid_wu_s") %>% mutate(antigen = "Spike model")
    data_i_all_rbd <- get_data_ind("sid_wu_rbd") %>% mutate(antigen = "RBD model")
    data_i_all <- bind_rows(data_i_all_s, data_i_all_rbd) %>%
        mutate(antigen = factor(antigen, levels = c("Spike model", "RBD model"))) 

    df_uncert_pop <- df_uncert_ind %>% filter(
                (state_variable == "Memory B-cell conc." & time < 80) |
                (state_variable == "Plasmablast conc." & time < 80) |
                (state_variable == "sVNT to ancestral variant" & time < 240)
                ) %>%
                group_by(time, state_variable, antigen) %>%
        mutate(antigen = factor(antigen, levels = c("Spike model", "RBD model"))) #%>%
             #   summarise(value = calculate_combined_mean(value, n), .lower = calculate_combined_ci_lower(value, .lower, .upper, n), .upper = calculate_combined_ci_upper(value, .lower, .upper, n))

    # Plot the population-level trajectories
    df_uncert_pop %>%
        ggplot() +
        stat_lineribbon(aes(x = time, y = value, fill = antigen), .width = c(0.5, 0.95),   alpha = 0.3) + 
     #   geom_line(aes(x = time, y = value, color = state_variable)) +
        geom_point(data = data_i_all %>% rename(state_variable = state), 
            aes(x = t, y = value, fill = antigen, shape = antigen), size = 3, alpha = 0.5) + theme_bw() +    
        scale_shape_manual(values = c(21, 22)) +
        scale_fill_manual(values = c(color_study_A, color_study_B)) +
        scale_color_manual(values = c(color_study_A, color_study_B)) +
        facet_grid(cols = vars(state_variable), rows = vars(antigen), scales = "free") + 
        labs(x = "Time after second dose (days)", y = "Conc. (%)/sVNT titre (log)", fill = "Data", color = "Posterior predictive interval") +
        guides(color = "none", fill = "none", shape = "none") + 
        theme(legend.position = "top") + 
        david_theme_B() +
        ggtitle(paste0("Model predictions to validation (unseen) dataset"))
    ggsave(here::here("outputs", "figs", "valid", "fig1C.png"), height = 8, width = 12)
}


#' Plot the population-level trajectories with the CRPS scores
plot_val_crps_ind <- function() {

   df_uncert_ind_s <- readRDS(here::here("outputs", "stanfit",  paste0("traj_val_",  "s", ".RDS"))) %>% mutate(antigen = "Spike model")
    df_uncert_ind_rbd <- readRDS(here::here("outputs", "stanfit",  paste0("traj_val_",  "rbd", ".RDS"))) %>% mutate(antigen = "RBD model")
    df_uncert_ind <- bind_rows(df_uncert_ind_s, df_uncert_ind_rbd) %>% filter(state %in% c("B", "P", "A")) %>%
        mutate(state_variable = recode(state, B = "Memory B-cell conc.", P = "Plasmablast conc.", A = "sVNT to ancestral variant"))

    # Get the calibration data
    data_i_all_s <- get_data_ind("sid_wu_s") %>% mutate(antigen = "Spike model")
    data_i_all_rbd <- get_data_ind("sid_wu_rbd") %>% mutate(antigen = "RBD model")
    data_i_all <- bind_rows(data_i_all_s, data_i_all_rbd)

    df_uncert_pop <- df_uncert_ind %>% filter(
                (state_variable == "Memory B-cell conc." & time < 80) |
                (state_variable == "Plasmablast conc." & time < 80) |
                (state_variable == "sVNT to ancestral variant" & time < 240)
                ) %>%
                group_by(time, state_variable, antigen) 


    state_names <-  c("Memory B-cell conc.", "Plasmablast conc.", "sVNT to ancestral variant")

    fit_val_compare <- data_i_all %>% rename(data_value = value) %>% left_join(df_uncert_pop %>% select(!state) %>% rename(t = time, state = state_variable))
    N <- fit_val_compare$PID %>% unique %>% length
     df_crps_all <- map_df(as.character(state_names),
        function(state_str) {
            #  "Memory B-cell conc."               "Plasmablast conc." "sVNT to ancestral variant (log)" 
           # bcell_mat <- y_hat_dist_data %>% ungroup %>% filter(!is.na(i), !is.na(data_value)) %>% select(.draw, PID, i, antigen, day_data, state, value, data_value) %>%
           #  state_str <- "sVNT to ancestral variant (log)"
            fit_val_compare_state <- fit_val_compare %>% filter(state == state_str) 
                
            mat_model_post <- fit_val_compare_state %>% select(antigen, PID, t, value) %>% pull(value) %>% matrix(nrow = 2 * N, byrow = TRUE) %>% t
            df_model_data <- fit_val_compare_state %>% select(antigen, PID, t, data_value) %>% unique 
            vec_model_data <- df_model_data %>% pull(data_value)

            crps_score_raw <- crps(x = mat_model_post, x2 = mat_model_post, y = vec_model_data) 
            crps_score <- 1 - (crps_score_raw$pointwise %>% exp)

            df_model_data$crps_score <- crps_score
            df_model_post_summary <- fit_val_compare_state %>% select(PID, t, value) %>% group_by(PID, t) %>% mean_qi(value)

            left_join(df_model_post_summary, df_model_data) %>% mutate(state = state_str)
        }
    ) 

    df_crps_all %>%
        mutate(antigen = factor(antigen, levels = c("Spike model", "RBD model"))) %>% 
        ggplot() + 
            geom_boxplot(aes(x = state, y = crps_score, fill = antigen), position = position_dodge(0.8), alpha = 0.3) + 
            theme_bw() +
            scale_fill_manual(values = c(color_study_A, color_study_B)) +
            david_theme_B() + 
            theme(legend.position = "top") + 
            labs(x = "", y = "Continuous Ranked Probability Score\n (CRPS)", fill = "", color = "") +        
             theme(strip.text.x = element_text(size = 24), axis.text = element_text(size = 20), axis.title.y = element_text(size = 22), 
              title = element_text(size = 20), legend.text = element_text(size = 20))   
    ggsave(here::here("outputs", "figs", "valid", "fig1D.png"), height = 10, width = 15)


    p1 <- df_crps_all %>% filter(antigen == "Spike model") %>% deanonymize_pid %>%
            ggplot() + 
                geom_linerange(aes(PID, y = value, ymin = .lower, ymax = .upper), size = 1) + 
                geom_point(aes(PID, y = value), color = "black", size = 3)  +
                geom_point(aes(PID, y = data_value, fill = crps_score), alpha = 0.7, shape = 22, color = "white", size = 3)  +
                facet_nested(vars(state), scales = "free") + theme_bw() +
                labs(x = "Person ID", y = "Value", fill = "Continous Ranked Probaility Score (CRPS)", color = "") + 
                scale_color_manual(labels = c( "Model fit"), values = c("black") ) + 
                # have low values red and high values gray in gradient
                scale_fill_continuous( low = "red", high = "gray")+ 
                theme(legend.position = "bottom") + 
                ggtitle(paste0("Model predictions for validation (unseen) " , "Ancestral spike",  " antigen") ) + 
                david_theme_A()
    p2 <- df_crps_all %>% filter(antigen == "RBD model") %>% deanonymize_pid %>%
            ggplot() + 
                geom_linerange(aes(PID, y = value, ymin = .lower, ymax = .upper), size = 1) + 
                geom_point(aes(PID, y = value), color = "black", size = 3)  +
                geom_point(aes(PID, y = data_value, fill = crps_score), alpha = 0.7, shape = 22, color = "white", size = 3)  +
                facet_nested(vars(state), scales = "free") + theme_bw() +
                labs(x = "Person ID", y = "Value", fill = "Continous Ranked Probaility Score (CRPS)", color = "") + 
                scale_color_manual(labels = c( "Model fit"), values = c("black") ) + 
                # have low values red and high values gray in gradient
                scale_fill_continuous( low = "red", high = "gray")+ 
                theme(legend.position = "bottom") + 
                ggtitle(paste0("Model predictions for validation (unseen) " , "Ancestral RBD",  " antigen") ) + 
                david_theme_A()
    p1 / p2
    ggsave(here::here("outputs", "figs", "valid", "figS2.png"), height = 15, width = 15)
}


#' Plot the population-level trajectories with the MRSE scores
plot_val_rmse_ind <- function() {

   df_uncert_ind_s <- readRDS(here::here("outputs", "stanfit",  paste0("traj_val_",  "s", ".RDS"))) %>% mutate(antigen = "Spike model")
    df_uncert_ind_rbd <- readRDS(here::here("outputs", "stanfit",  paste0("traj_val_",  "rbd", ".RDS"))) %>% mutate(antigen = "RBD model")
    df_uncert_ind <- bind_rows(df_uncert_ind_s, df_uncert_ind_rbd) %>% filter(state %in% c("B", "P", "A")) %>%
        mutate(state_variable = recode(state, B = "Memory B-cell conc.", P = "Plasmablast conc.", A = "sVNT to ancestral variant"))

    # Get the calibration data
    data_i_all_s <- get_data_ind("sid_wu_s") %>% mutate(antigen = "Spike model")
    data_i_all_rbd <- get_data_ind("sid_wu_rbd") %>% mutate(antigen = "RBD model")
    data_i_all <- bind_rows(data_i_all_s, data_i_all_rbd)

    df_uncert_pop <- df_uncert_ind %>% filter(
                (state_variable == "Memory B-cell conc." & time < 80) |
                (state_variable == "Plasmablast conc." & time < 80) |
                (state_variable == "sVNT to ancestral variant" & time < 240)
                ) %>%
                group_by(time, state_variable, antigen) 


    state_names <-  c("Memory B-cell conc.", "Plasmablast conc.", "sVNT to ancestral variant")

    fit_val_compare <- data_i_all %>% rename(data_value = value) %>% left_join(df_uncert_pop %>% select(!state) %>% rename(t = time, state = state_variable))
    N <- fit_val_compare$PID %>% unique %>% length
     df_rmse_all <- map_df(as.character(state_names),
        function(state_str) {
                        #  "Memory B-cell conc."               "Plasmablast conc." "sVNT to ancestral variant (log)" 
           # bcell_mat <- y_hat_dist_data %>% ungroup %>% filter(!is.na(i), !is.na(data_value)) %>% select(.draw, PID, i, antigen, day_data, state, value, data_value) %>%
           #  state_str <- "sVNT to ancestral variant (log)"
            fit_val_compare_state <- fit_val_compare %>% filter(state == state_str) 
                
            mat_model_post <- fit_val_compare_state %>% select(antigen, PID, t, value) %>% pull(value) %>% matrix(nrow = 2 * N, byrow = TRUE) %>% t
            df_model_data <- fit_val_compare_state %>% select(antigen, PID, t, data_value) %>% unique 
            vec_model_data <- df_model_data %>% pull(data_value)

            pred_mean <- colMeans(mat_model_post)
            rmse_score <- sqrt((vec_model_data - pred_mean)^2)

            df_model_data$rmse_score <- rmse_score
        
            df_model_post_summary <- fit_val_compare_state %>% select(PID, t, value) %>% group_by(PID, t) %>% mean_qi(value)

            left_join(df_model_post_summary, df_model_data) %>% mutate(state = state_str)
        }
    ) 

    df_rmse_all %>%
        mutate(antigen = factor(antigen, levels = c("Spike model", "RBD model"))) %>% 
        ggplot() + 
            geom_boxplot(aes(x = state, y = rmse_score, fill = antigen), position = position_dodge(0.8), alpha = 0.3) + 
            theme_bw() +
            scale_fill_manual(values = c(color_study_A, color_study_B)) +
            david_theme_B() + 
            theme(legend.position = "top") + 
            labs(x = "", y = "Root Mean Squared Error\n (RMSE)", fill = "", color = "") +        
             theme(strip.text.x = element_text(size = 24), axis.text = element_text(size = 20), axis.title.y = element_text(size = 22), 
              title = element_text(size = 20), legend.text = element_text(size = 20))   
    ggsave(here::here("outputs", "figs", "valid", "fig1D_RMSE.png"), height = 10, width = 15)


    p1 <- df_rmse_all %>% filter(antigen == "Spike model") %>% deanonymize_pid %>%
            ggplot() + 
                geom_linerange(aes(PID, y = value, ymin = .lower, ymax = .upper), size = 1) + 
                geom_point(aes(PID, y = value), color = "black", size = 3)  +
                geom_point(aes(PID, y = data_value, fill = rmse_score), alpha = 0.7, shape = 22, color = "white", size = 3)  +
                facet_nested(vars(state), scales = "free") + theme_bw() +
                labs(x = "Person ID", y = "Value", fill = "Root Mean Squared Error (RMSE)", color = "") + 
                scale_color_manual(labels = c( "Model fit"), values = c("black") ) + 
                # have low values red and high values gray in gradient
                scale_fill_continuous( low = "red", high = "gray")+ 
                theme(legend.position = "bottom") + 
                ggtitle(paste0("Model predictions for validation (unseen) " , "Ancestral spike",  " antigen") ) + 
                david_theme_A()
    p2 <- df_rmse_all %>% filter(antigen == "RBD model") %>% deanonymize_pid %>%
            ggplot() + 
                geom_linerange(aes(PID, y = value, ymin = .lower, ymax = .upper), size = 1) + 
                geom_point(aes(PID, y = value), color = "black", size = 3)  +
                geom_point(aes(PID, y = data_value, fill = rmse_score), alpha = 0.7, shape = 22, color = "white", size = 3)  +
                facet_nested(vars(state), scales = "free") + theme_bw() +
                labs(x = "Person ID", y = "Value", fill = "Root Mean Squared Error (RMSE)", color = "") + 
                scale_color_manual(labels = c( "Model fit"), values = c("black") ) + 
                # have low values red and high values gray in gradient
                scale_fill_continuous( low = "red", high = "gray")+ 
                theme(legend.position = "bottom") + 
                ggtitle(paste0("Model predictions for validation (unseen) " , "Ancestral RBD",  " antigen") ) + 
                david_theme_A()
    p1 / p2
    ggsave(here::here("outputs", "figs", "valid", "figS2_RMSE.png"), height = 15, width = 15)
}


#' Plot the individual-level trajectories for the validations data given the stanfit object for each individual
plot_val_traj_ind <- function() {

   # Get the trajectories for each individual
    df_uncert_ind_s <- readRDS(here::here("outputs", "stanfit",  paste0("traj_val_",  "s", ".RDS")))  %>% mean_qi(value) %>% mutate(antigen = "Ancestral spike")
    df_uncert_ind_rbd <- readRDS(here::here("outputs", "stanfit",  paste0("traj_val_",  "rbd", ".RDS")))  %>% mean_qi(value) %>% mutate(antigen = "Ancestral RBD")
    df_uncert_ind <- bind_rows(df_uncert_ind_s, df_uncert_ind_rbd) %>% filter(state %in% c("B", "P", "A")) %>%
        mutate(state_variable = recode(state, B = "Memory B-cell conc.", P = "Plasmablast conc.", A = "sVNT to ancestral variant")) %>% mutate(n = 4000) %>%
        mutate(antigen = factor(antigen, levels = c("Ancestral spike", "Ancestral RBD"))) 

    # Get the calibration data
    data_i_all_s <- get_data_ind("sid_wu_s") %>% mutate(antigen = "Ancestral spike")
    data_i_all_rbd <- get_data_ind("sid_wu_rbd") %>% mutate(antigen = "Ancestral RBD")
    data_i_all <- bind_rows(data_i_all_s, data_i_all_rbd) %>%
        mutate(antigen = factor(antigen, levels = c("Ancestral spike", "Ancestral RBD"))) 

    p1 <- df_uncert_ind %>% filter(antigen == "Ancestral spike", time < 220) %>% 
        deanonymize_pid %>%
        ggplot() +
        geom_ribbon(aes(x = time, y = value, ymin = .lower, ymax = .upper, fill = state_variable), 
        alpha = 0.3) + 
        geom_line(aes(x = time, y = value, color = state_variable)) + 
        geom_point(data = data_i_all  %>% filter(antigen == "Ancestral spike") %>% deanonymize_pid, aes(x = t, y = value, fill = state), shape = 21, size = 3, alpha = 0.5) + theme_bw() +    
        xlim(0, 80) +
        facet_wrap(vars(PID)) + labs(x = "Days post-vaccination", y = "Conc./sVNT titre", fill = "Data", color = "Posterior predictive interval") + 
        theme(legend.position = "top") + 
        ggtitle(paste0("Individual-level immune trajectories in response to ",  "Ancestral spike", " antigen on validation data") )

    p2 <- df_uncert_ind %>% filter(antigen == "Ancestral RBD", time < 220) %>% deanonymize_pid %>%
        ggplot() +
        geom_ribbon(aes(x = time, y = value, ymin = .lower, ymax = .upper, fill = state_variable), 
        alpha = 0.3) + 
        geom_line(aes(x = time, y = value, color = state_variable)) + 
        geom_point(data = data_i_all %>% filter(antigen == "Ancestral RBD") %>% deanonymize_pid, aes(x = t, y = value, fill = state), shape = 21, size = 3, alpha = 0.5) + theme_bw() +    
        xlim(0, 80) +
        facet_wrap(vars(PID)) + labs(x = "Days post-vaccination", y = "Conc./sVNT titre", fill = "Data", color = "Posterior predictive interval") + 
        theme(legend.position = "top") + 
        ggtitle(paste0("Individual-level immune trajectories in response to ",  "Ancestral RBD", " antigen on validation data") )
   # ggsave(here::here("outputs", "figs", "manu_A", foldername, "fits", "fig1_ind_traj_XX.png"), width = 15, height = 12, units = "in")
    #p1 <- p4A  + plot_layout(heights = c(1, 2)) + plot_annotation(tag_levels = "A") + plot_layout(guides = "collect") &  theme(legend.position = "top") 
    p1 / p2 + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") & theme(legend.position='top') 
    ggsave(here::here("outputs", "figs", "valid", "figS4.png"), height = 20)
}


#' This is old code but maybe useful
compare_cprs_antibodies <- function(stanfit, fig_folder) {   
  
    if (fig_folder == "sid_wu_s") {
        antigen_str <- "Ancestral spike"
        antigen_file <- "s"
    } else {
        antigen_str <- "Ancestral RBD"
        antigen_file <- "rbd"
    }
   # post_sample_cross <- readRDS(file = here::here("outputs", "stanfit", "manu_A", paste0(stanfit_name, ".RDS")) )
    data_meta <- readRDS(here::here("outputs", "data_clean", fig_folder, "df_meta.RDS")) %>% 
        mutate(time_d1_i = recode(time_d1, "<28 days" = "1", "28+ days" = "2")) %>% 
        mutate(vax_type_i = recode(vax_type, "BNT162b2" = "1", "ChAdOx1" = "2"))


    data_model <- readRDS(here::here("outputs", "data_clean", fig_folder, "df_model_data.RDS"))
    data_listA <- readRDS(here::here("outputs", "data_clean", fig_folder, "list_stan_data.RDS"))
    data_listA_cal <- readRDS(here::here("outputs", "data_clean", "nih_wu_s", "list_stan_data.RDS"))

    data_listA_cal$x_vac %>% table

    PIDs <- data_meta$PID
    # theta_t and theta_a
    df_compare_time_age <- map_df(list("theta_t", "theta_a"), 
        function(fitted_theta) {
            param_expr = expr((!!sym(fitted_theta))[j, i, s])
            param_uncert <- stanfit$draws(
                c(fitted_theta,   
                    "prop_asc", "deltaAg", "deltaP", "deltaL", "deltaA", "deltaG")) %>%
                        spread_draws(!!param_expr, prop_asc, deltaAg, deltaP, deltaL, deltaA, deltaG) %>%
                mutate(deltaM = 0.001, .before = "deltaP") %>%
                mutate(deltaAg = 1 / (deltaAg * 30), deltaP = 1 / deltaP, deltaL = 1 / deltaL, deltaA = 1 / deltaA, deltaG = 1 / deltaG) %>% 
                pivot_wider(names_from = "s", values_from = fitted_theta) %>% 
                select(j:.draw, `1`:`4`, prop_asc, deltaAg:deltaG) %>% 
                rename(t = j, v = i) %>% rename(param1 = `1`, param2 = `2`, param3 = `3`, param4 = `4`)

            param_uncert_sum <- param_uncert %>% ungroup %>% summarise(across(param1:deltaG, ~mean(.x)), .by = c("t", "v")) 

            M <- 100
            df_uncert_all <- map_df(1:length(PIDs), 
                function(s) {
                    init_i <- data_listA$y0_init[s, ]
                    names(init_i) <- c("B", "P", "A1", "G", "L", "A2")

                    param_uncert_sum_i <- param_uncert %>% filter(t == data_meta[s, ]$time_d1_i, v == data_meta[s, ]$vax_type_i)

                    df_uncert <- map_df(1:M,
                        ~ode_results_vac(param_uncert_sum_i[.x, 6:16], 80, init_i) %>% as.data.frame  )  %>%
                        select(time, B, P, A1, A2) %>% mutate(A = A1 + A2) %>% 
                        mutate(i = map(1:M, ~rep(.x, 81)) %>% unlist) %>% 
                        mutate(PID = PIDs[s])
                }
            )
            df_uncert_all_add <- df_uncert_all %>% left_join(data_meta)

            df_data_A <- data_model %>% select(PID, time_post_exp, IC50hs, vax_type, time_d1, age_cat) %>% filter(time_post_exp != 0)
            df_crps <- map_df(1:length(df_data_A$PID),
                ~df_uncert_all_add %>% filter(PID == df_data_A[[.x, 1]], time == df_data_A[[.x, 2]]) %>% mutate(data = df_data_A[[.x, 3]])
            )
            
            df_cprs_data <- df_crps %>% select(PID, data) %>% unique %>% pull(data)
            df_cprs_mat <- df_crps %>% pull(A) %>% matrix(c(length(df_data_A$PID)), byrow = TRUE) %>% t

            library(loo)
            crps_score <- crps(x = df_cprs_mat, x2 = df_cprs_mat, y = df_cprs_data) 

            data.frame(
                PID = df_data_A$PID,
                crps = crps_score$pointwise,
                mean =  crps_score$estimate[[1]],
                sd =  crps_score$estimate[[2]],
                type = fitted_theta
            )

           # crps_score$estimates %>% as.data.frame %>% t %>% as.data.frame %>% mutate(type = !!fitted_theta, alt_score = exp(crps_score$pointwise) %>% mean)
    
        }
    )

    param_uncert <- stanfit$draws(
        c("theta_v",   
            "prop_asc", "deltaAg", "deltaP", "deltaL", "deltaA", "deltaG")) %>%
                spread_draws(theta_v[i, s], prop_asc, deltaAg, deltaP, deltaL, deltaA, deltaG) %>%
        mutate(deltaM = 0.001, .before = "deltaP") %>%
        mutate(deltaAg = 1 / (deltaAg * 30), deltaP = 1 / deltaP, deltaL = 1 / deltaL, deltaA = 1 / deltaA, deltaG = 1 / deltaG) %>% 
        pivot_wider(names_from = "s", values_from = "theta_v") %>% 
        select(i:.draw, `1`:`4`, prop_asc, deltaAg:deltaG) %>% 
        rename(v = i) %>% rename(param1 = `1`, param2 = `2`, param3 = `3`, param4 = `4`)

    param_uncert_sum <- param_uncert %>% ungroup %>% summarise(across(param1:deltaG, ~mean(.x)), .by = c("v")) 

    M <- 100
    df_uncert_all <- map_df(1:length(PIDs), 
        function(s) {
            init_i <- data_listA$y0_init[s, ]
            names(init_i) <- c("B", "P", "A1", "G", "L", "A2")

            param_uncert_sum_i <- param_uncert %>% filter( v == data_meta[s, ]$vax_type_i)

            df_uncert <- map_df(1:M,
                ~ode_results_vac(param_uncert_sum_i[.x, 5:15], 80, init_i) %>% as.data.frame  )  %>%
                select(time, B, P, A1, A2) %>% mutate(A = A1 + A2) %>% 
                mutate(i = map(1:M, ~rep(.x, 81)) %>% unlist) %>% 
                mutate(PID = PIDs[s])
        }
    )
    df_uncert_all_add <- df_uncert_all %>% left_join(data_meta)

    df_data_A <- data_model %>% select(PID, time_post_exp, IC50hs, vax_type, time_d1, age_cat) %>% filter(time_post_exp != 0)
    df_crps <- map_df(1:length(df_data_A$PID),
        ~df_uncert_all_add %>% filter(PID == df_data_A[[.x, 1]], time == df_data_A[[.x, 2]]) %>% mutate(data = df_data_A[[.x, 3]])
    )
    
    df_cprs_data <- df_crps %>% select(PID, data) %>% unique %>% pull(data)
    df_cprs_mat <- df_crps %>% pull(A) %>% matrix(c(length(df_data_A$PID)), byrow = TRUE) %>% t
    crps_score <- crps(x = df_cprs_mat, x2 = df_cprs_mat, y = df_cprs_data) 


    df_compare_vac_time_age <- bind_rows(
        df_compare_time_age,
        data.frame(
            PID = df_data_A$PID,
            crps = crps_score$pointwise,
            mean =  crps_score$estimate[[1]],
            sd =  crps_score$estimate[[2]],
            type = "theta_v"
        )    
    )


    param_uncert <- stanfit$draws(
        c("theta_main",   
            "prop_asc", "deltaAg", "deltaP", "deltaL", "deltaA", "deltaG")) %>%
                spread_draws(theta_main[s], prop_asc, deltaAg, deltaP, deltaL, deltaA, deltaG) %>%
        mutate(deltaM = 0.001, .before = "deltaP") %>%
        mutate(deltaAg = 1 / (deltaAg * 30), deltaP = 1 / deltaP, deltaL = 1 / deltaL, deltaA = 1 / deltaA, deltaG = 1 / deltaG) %>% 
        pivot_wider(names_from = "s", values_from = "theta_v") %>% 
        select(i:.draw, `1`:`4`, prop_asc, deltaAg:deltaG) %>% 
        rename(v = i) %>% rename(param1 = `1`, param2 = `2`, param3 = `3`, param4 = `4`)

    param_uncert_sum <- param_uncert %>% ungroup %>% summarise(across(param1:deltaG, ~mean(.x)), .by = c("v")) 


   param_uncert_sum_none <- param_uncert %>% ungroup %>% summarise(across(param1:deltaG, ~weighted.mean(.x, c( 0.8048, 1 -  0.8048))), .by = .draw)


    M <- 100
    df_uncert_all_none <- map_df(1:length(PIDs), 
        function(s) {
            init_i <- data_listA$y0_init[s, ]
            names(init_i) <- c("B", "P", "A1", "G", "L", "A2")


            df_uncert <- map_df(1:M,
                ~ode_results_vac(param_uncert_sum_none[.x, 2:12], 80, init_i) %>% as.data.frame  )  %>%
                select(time, B, P, A1, A2) %>% mutate(A = A1 + A2) %>% 
                mutate(i = map(1:M, ~rep(.x, 81)) %>% unlist) %>% 
                mutate(PID = PIDs[s])
        }
    )
    df_uncert_all_add_none <- df_uncert_all_none %>% left_join(data_meta)

    df_data_A <- data_model %>% select(PID, time_post_exp, IC50hs, vax_type, time_d1, age_cat) %>% filter(time_post_exp != 0)
    df_crps_none <- map_df(1:length(df_data_A$PID),
        ~df_uncert_all_add_none %>% filter(PID == df_data_A[[.x, 1]], time == df_data_A[[.x, 2]]) %>% mutate(data = df_data_A[[.x, 3]])
    )

    df_cprs_data <- df_crps_none %>% select(PID, data) %>% unique %>% pull(data)
    df_cprs_mat <- df_crps_none %>% pull(A) %>% matrix(c(length(df_data_A$PID)), byrow = TRUE) %>% t
    crps_score_none <- crps(x = df_cprs_mat, x2 = df_cprs_mat, y = df_cprs_data) 

    df_compare_vac_time_age_all <- bind_rows(
        df_compare_vac_time_age,
        data.frame(
                PID = df_data_A$PID,
                crps = crps_score_none$pointwise,
                mean =  crps_score_none$estimate[[1]],
                sd =  crps_score_none$estimate[[2]],
                type = "theta"
            )    
    )

    df_compare_vac_time_age_all %>% select( mean, sd, type) %>% unique %>% 
        mutate(lb = mean - 2 * sd, ub = mean + 2 * sd) %>%
        mutate(desc = c("Vaccine type and dose timing", "Vaccine type and age", "Vaccine type", "None")) %>%
        mutate(desc = factor(desc, levels = c("Vaccine type and dose timing", "Vaccine type and age", "Vaccine type", "None"))) %>%
        ggplot() + 
            guides(color = "none") +
            geom_linerange(aes(y = desc, xmin = lb,  xmax = ub, color = type), size = 1) +
            geom_point(aes(y = desc, x = mean, color = type), size = 5) + theme_bw()  + 
            ggtitle("CRPS for different host factors in mechanistic model") + 
            labs(x = "Average CRPS across all individuals", y = "")


    ggsave(here::here("outputs", "figs", "valid", paste0("CRPS_score_",  antigen_file, ".png")), width = 15, height = 12, units = "in")

}