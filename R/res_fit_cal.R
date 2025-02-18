
#' Calculate the trajectories for the calibration given the stanfit object for each individual
#' @param stanfit The stanfit object
#' @param foldername The foldername
cal_traj_ind <- function(stanfit, foldername) {
    
    data_meta <- readRDS(here::here("outputs", "data_clean", foldername, "df_meta.RDS"))
    data_listA <- readRDS(here::here("outputs", "data_clean", foldername, "list_stan_data.RDS"))

    PIDs <- data_meta$PID 

    param_uncert <- stanfit$draws(
        c("param1", "v_1", "sigmav1", "t_1", "sigmat1", "a_1", "sigmaa1", "ind_1", "sigmai1", 
        "param2",
        "param3", "v_3", "sigmav3", "t_3", "sigmat3", "a_3", "sigmaa3", "ind_3", "sigmai3", 
        "param4", "v_4", "sigmav4", "t_4", "sigmat4", "a_4", "sigmaa4", "ind_4", "sigmai4",   
            "prop_asc", "deltaAg", "deltaP", "deltaL", "deltaA", "deltaG")) %>%
        spread_draws(
            param1, v_1[v], sigmav1, t_1[t], sigmat1, a_1[a], sigmaa1, ind_1[i], sigmai1, 
            param2,
            param3, v_3[v], sigmav3, t_3[t], sigmat3, a_3[a], sigmaa3, ind_3[i], sigmai3,
            param4, v_4[v], sigmav4, t_4[t], sigmat4, a_4[a], sigmaa4, ind_4[i], sigmai4,
            prop_asc, deltaAg, deltaP, deltaL, deltaA, deltaG) %>% 
        mutate(
            param1 = inv.logit(param1 + v_1*sigmav1 + t_1*sigmat1 + a_1*sigmaa1 + ind_1*sigmai1 ) * 2,
            param2 = inv.logit(param2),
            param3 = inv.logit(param3 + v_3*sigmav3 + t_3*sigmat3 + a_3*sigmaa3 + ind_3*sigmai3 ) * 6,
            param4 = inv.logit(param4 + v_4*sigmav4 + t_4*sigmat4 + a_4*sigmaa4 + ind_4*sigmai4 ) * 0.5, .keep = "unused") %>%
        mutate(deltaM = 0.001, .before = "deltaP") %>%
        mutate(deltaAg = 1 / (deltaAg * 30), deltaP = 1 / deltaP, deltaL = 1 / deltaL, deltaA = 1 / deltaA, deltaG = 1 / deltaG)

    param_uncert <- param_uncert %>% select(.iteration, v, t, a, i, param1, param2:deltaG)
    param_uncert_sum <- param_uncert %>% ungroup %>% summarise(across(param1:deltaG, ~mean(.x)), .by = c("v", "t", "i", "a")) 

    pid_names <- PIDs %>% unique
    N <- length(pid_names)

    M <- 100

    df_uncert_all <- map_df( 1:N,
    function(ind) {

        init_i <- data_listA$y0_init[ind, ]#https://file+.vscode-resource.vscode-cdn.net/Users/davidhodgson/Dropbox/Mac%20%283%29/Documents/research/nih/covid_bcell_published/outputs/figs/post/fig1_fit_ind_A.png?version%3D1717687880997
        names(init_i) <- c("B", "P", "A1", "G", "L", "A2")


        param_uncert_sum_i <- param_uncert %>% filter(i == ind, a == data_listA$x_age[ind], t == data_listA$x_time[ind], v == data_listA$x_vac[ind])

        df_uncert <- map_df(1:M,
            ~ode_results_vac(param_uncert_sum_i[.x, 6:16], 365, init_i) %>% as.data.frame  )  %>%
            select(time, B, P, A1, A2) %>% mutate(A = A1 + A2, .keep = "unused") %>% 
        pivot_longer(!time, names_to = "state_variable", values_to = "value") %>%
            group_by(time, state_variable) %>% mean_qi(value) %>% mutate(PID = PIDs[ind], n = M)
        }
    )

    relabel_state_variable <- c("A" = "sVNT to ancestral variant", "P" = "Plasmablast conc.", "B" = "Memory B-cell conc.")
    df_uncert_all_labs <- df_uncert_all %>% mutate(state_variable = recode(state_variable, !!!relabel_state_variable) )
    df_uncert_all_labs %>% left_join(data_model %>% select(PID, vax_type, time_d1,  age_cat   ) )
    df_uncert_all_labs
}


#' Run the full posterior trajectories
#' @param stanfit The stanfit object
#' @param foldername The foldername
run_full_posterior_tracjectories <- function(stanfit, foldername) {
    if (foldername == "nih_wu_s") {
        antigen_str <- "Ancestral spike"
        antigen_file <- "s"
    } else {
        antigen_str <- "Ancestral RBD"
        antigen_file <- "rbd"
    }
    df_uncert_all_labs_s <- cal_traj_ind(stanfit, foldername)
    saveRDS(df_uncert_all_labs_s, here::here("outputs", "stanfit",  paste0("traj_ind_",  antigen_file, ".RDS")))
}

#' Plot the population-level trajectories
plot_cal_traj_pop <- function() {

    # Get the trajectories for each individual
    df_uncert_ind_s <- readRDS(here::here("outputs", "stanfit",  paste0("traj_ind_",  "s", ".RDS"))) %>% mutate(antigen = "Spike model")
    df_uncert_ind_rbd <- readRDS(here::here("outputs", "stanfit",  paste0("traj_ind_",  "rbd", ".RDS"))) %>% mutate(antigen = "RBD model")
    df_uncert_ind <- bind_rows(df_uncert_ind_s, df_uncert_ind_rbd)

    # Get the calibration data
    data_i_all_s <- get_data_ind("nih_wu_s") %>% mutate(antigen = "Spike model")
    data_i_all_rbd <- get_data_ind("nih_wu_rbd") %>% mutate(antigen = "RBD model")
    data_i_all <- bind_rows(data_i_all_s, data_i_all_rbd) %>%
        mutate(antigen = factor(antigen, levels = c("Spike model", "RBD model"))) 

    # Combine the individual-level trajectories to population-level trajectories
    df_uncert_pop <- df_uncert_ind %>% filter(
            (state_variable == "Memory B-cell conc." & time < 80) |
            (state_variable == "Plasmablast conc." & time < 80) |
            (state_variable == "sVNT to ancestral variant" & time < 240)
            ) %>%
            group_by(time, state_variable, antigen) %>%
            mutate(antigen = factor(antigen, levels = c("Spike model", "RBD model"))) 
           # summarise(value = calculate_combined_mean(value, n), .lower = calculate_combined_ci_lower(value, .lower, .upper, n), .upper = calculate_combined_ci_upper(value, .lower, .upper, n))

    # Plot the population-level trajectories
    df_uncert_pop %>%
        ggplot() +
        stat_lineribbon(aes(x = time, y = value, fill = antigen), .width = c(0.5, 0.95),   alpha = 0.3) + 
      #  geom_line(aes(x = time, y = value, color = state_variable)) +
        geom_point(data = data_i_all %>% rename(state_variable = state), 
            aes(x = t, y = value, fill = antigen, shape = antigen), size = 3, alpha = 0.5) + theme_bw() +    
        scale_shape_manual(values = c(21, 22)) +
        scale_fill_manual(values = c(color_study_A, color_study_B)) +
        scale_color_manual(values = c(color_study_A, color_study_B)) +
        facet_grid(cols = vars(state_variable), rows = vars(antigen), scales = "free") + 
        labs(x = "Time after second dose (days)", y = "Conc. (%)/sVNT (log)", fill = "Data", color = "Posterior predictive interval") +
        guides(color = "none", fill = "none", shape = "none") + 
        theme(legend.position = "top") + 
        david_theme_B() +
        ggtitle(paste0("Model fits to calibration dataset"))
    ggsave(here::here("outputs", "figs", "post_pred", "fig1A.png"), height = 8, width = 12)
}

#' Plot the population-level trajectories with the CRPS calculated
plot_cal_crps_ind <- function() {

    # Get the calibration data
    data_i_all_s <- get_data_ind("nih_wu_s") %>% mutate(antigen = "Spike model")
    data_i_all_rbd <- get_data_ind("nih_wu_rbd") %>% mutate(antigen = "RBD model")
    data_i_all <- bind_rows(data_i_all_s, data_i_all_rbd) %>% deanonymize_pid() %>% 
        group_by(antigen, state, PID) %>% mutate(t_id = row_number()) %>%
        rename(data_value = value)

    state_names <-  c("Memory B-cell conc.", "Plasmablast conc.", "sVNT to ancestral variant")
    state_id_data <- c("y_hat_1", "y_hat_2", "y_hat_3")
    names(state_names) <- state_id_data

    stanfit_s <- readRDS(file = here::here("outputs", "stanfit", paste0("nih_vac_wu_s_i", ".RDS")) )
    stanfit_rbd <- readRDS(file = here::here("outputs", "stanfit", paste0("nih_vac_wu_rbd_i", ".RDS")) )

    y_hat_dist_s <- stanfit_s$draws(c("y_hat")) %>% spread_draws(y_hat[n, t, s]) %>% rename(PID = n, t_id = t) %>%
        mutate(PID = factor(PID, levels = 1:41)) %>%  filter(s %in% c(1, 2, 3, 6)) %>%
            pivot_wider(names_from = "s", values_from = c("y_hat")) %>% 
            rename(y_hat_1 = `1`, y_hat_2 = `2`, y_hat_3 = `3`, y_hat_6 = `6`) %>%
            mutate( y_hat_3 = y_hat_3 + y_hat_6, .keep = "unused") %>% 
            pivot_longer(y_hat_1:y_hat_3, names_to = "state", values_to = "value") %>%
            mutate(state = recode(state, !!!state_names)) %>% mutate(antigen = "Spike model")

    y_hat_dist_rbd <- stanfit_rbd$draws(c("y_hat")) %>% spread_draws(y_hat[n, t, s]) %>% rename(PID = n, t_id = t) %>%
        mutate(PID = factor(PID, levels = 1:41)) %>%  filter(s %in% c(1, 2, 3, 6)) %>%
            pivot_wider(names_from = "s", values_from = c("y_hat")) %>% 
            rename(y_hat_1 = `1`, y_hat_2 = `2`, y_hat_3 = `3`, y_hat_6 = `6`) %>%
            mutate( y_hat_3 = y_hat_3 + y_hat_6, .keep = "unused") %>% 
            pivot_longer(y_hat_1:y_hat_3, names_to = "state", values_to = "value") %>%
            mutate(state = recode(state, !!!state_names)) %>% mutate(antigen = "RBD model")
    y_hat_dist <- bind_rows(y_hat_dist_s, y_hat_dist_rbd) 
    data_i_all %>% head

    y_hat_dist_data <- data_i_all %>% left_join(y_hat_dist %>% filter(!is.na(value))) %>% ungroup %>% 
        mutate(key = paste(PID, t_id, sep = "_")) %>%
        mutate(key = factor(key, levels = unique(key))) %>%
        mutate(key = as.numeric(key)) 

    df_crps_all <- map_df(as.character(state_names),
        function(state_str) {
            #  "Memory B-cell conc."               "Plasmablast conc." "sVNT to ancestral variant (log)" 
           # bcell_mat <- y_hat_dist_data %>% ungroup %>% filter(!is.na(i), !is.na(data_value)) %>% select(.draw, PID, i, antigen, day_data, state, value, data_value) %>%
           #  state_str <- "sVNT to ancestral variant (log)"
            y_hat_dist_data_state <- y_hat_dist_data %>% filter(state == state_str) 
                
            mat_model_post <- y_hat_dist_data_state %>% select(antigen, key, PID, t, value) %>% pull(value) %>% matrix(ncol = 4000, byrow = TRUE) %>% t
            df_model_data <- y_hat_dist_data_state %>% select(antigen, key, PID, t, data_value) %>% unique 
            vec_model_data <- df_model_data %>% pull(data_value)

            crps_score_raw <- crps(x = mat_model_post, x2 = mat_model_post, y = vec_model_data) 
            crps_score <- 1 - (crps_score_raw$pointwise %>% exp)

            df_model_data$crps_score <- crps_score
            df_model_post_summary <- y_hat_dist_data_state %>% select(key, PID, t, value) %>% group_by(PID, key, t) %>% mean_qi(value)

            left_join(df_model_post_summary, df_model_data) %>% mutate(state = state_str)
        }
    )  %>% mutate(VisitN = cut(t, breaks = c(-1, 13, 22, 300), labels = c("2–13 days", "14–21 days", "67–210 days") ) ) %>%
         mutate(VisitN = factor(VisitN, levels = c("2–13 days", "14–21 days", "67–210 days")) ) %>% filter(!is.na(VisitN)) %>%
        mutate(antigen = factor(antigen, levels = c("Spike model", "RBD model"))) 

    df_crps_all %>% 
        ggplot() + 
            geom_boxplot(aes(x = state, y = crps_score, fill = antigen), position = position_dodge(0.8), alpha = 0.3) + 
            theme_bw() +
            scale_fill_manual(values = c(color_study_A, color_study_B)) + 
            david_theme_B() + 
            theme(legend.position = "top") + 
            labs(x = "", y = "Continuous Ranked Probability Score\n (CRPS)", fill = "", color = "") +        
             theme(strip.text.x = element_text(size = 24), axis.text = element_text(size = 20), axis.title.y = element_text(size = 22), 
              title = element_text(size = 20), legend.text = element_text(size = 20))   
    ggsave(here::here("outputs", "figs", "post_pred", "fig1B.png"), height = 10, width = 15 )


    p1 <- df_crps_all %>% filter(antigen == "Spike model") %>%
            ggplot() + 
                geom_linerange(aes(PID, y = value, ymin = .lower, ymax = .upper), size = 1) + 
                geom_point(aes(PID, y = value), color = "black", size = 3)  +
                geom_point(aes(PID, y = data_value, fill = crps_score), alpha = 0.7, shape = 22, color = "white", size = 3)  +
                facet_nested(VisitN ~ state, scales = "free") + theme_bw() +
                labs(x = "Person ID", y = "Value", fill = "Continous Ranked Probaility Score (CRPS)", color = "") + 
                scale_color_manual(labels = c( "Model fit"), values = c("black") ) + 
                # have low values red and high values gray in gradient
                scale_fill_continuous( low = "red", high = "gray")+ 
                theme(legend.position = "bottom") + 
                ggtitle(paste0("Model fits to calibration data for " , "Ancestral spike",  " antigen") ) + 
                david_theme_A() + theme(axis.text.x = element_blank())
    p2 <- df_crps_all %>% filter(antigen == "RBD model") %>%
            ggplot() + 
                geom_linerange(aes(PID, y = value, ymin = .lower, ymax = .upper), size = 1) + 
                geom_point(aes(PID, y = value), color = "black", size = 3)  +
                geom_point(aes(PID, y = data_value, fill = crps_score), alpha = 0.7, shape = 22, color = "white", size = 3)  +
                facet_nested(VisitN ~ state, scales = "free") + theme_bw() +
                labs(x = "Person ID", y = "Value", fill = "Continous Ranked Probaility Score (CRPS)", color = "") + 
                scale_color_manual(labels = c( "Model fit"), values = c("black") ) + 
                # have low values red and high values gray in gradient
                scale_fill_continuous( low = "red", high = "gray")+ 
                theme(legend.position = "bottom") + 
                ggtitle(paste0("Model fits for calibration data" , " Ancestral RBD",  " antigen") ) + 
                david_theme_A() + theme(axis.text.x = element_blank())
    p1 / p2
    ggsave(here::here("outputs", "figs", "post_pred", "figS1.png"), height = 15, width = 15)

}


#' Plot the population-level trajectories with the CRPS calculated
plot_cal_rmse_ind <- function() {

    # Get the calibration data
    data_i_all_s <- get_data_ind("nih_wu_s") %>% mutate(antigen = "Spike model")
    data_i_all_rbd <- get_data_ind("nih_wu_rbd") %>% mutate(antigen = "RBD model")
    data_i_all <- bind_rows(data_i_all_s, data_i_all_rbd) %>% deanonymize_pid() %>% 
        group_by(antigen, state, PID) %>% mutate(t_id = row_number()) %>%
        rename(data_value = value)

    state_names <-  c("Memory B-cell conc.", "Plasmablast conc.", "sVNT to ancestral variant")
    state_id_data <- c("y_hat_1", "y_hat_2", "y_hat_3")
    names(state_names) <- state_id_data

    stanfit_s <- readRDS(file = here::here("outputs", "stanfit", paste0("nih_vac_wu_s_i", ".RDS")) )
    stanfit_rbd <- readRDS(file = here::here("outputs", "stanfit", paste0("nih_vac_wu_rbd_i", ".RDS")) )

    y_hat_dist_s <- stanfit_s$draws(c("y_hat")) %>% spread_draws(y_hat[n, t, s]) %>% rename(PID = n, t_id = t) %>%
        mutate(PID = factor(PID, levels = 1:41)) %>%  filter(s %in% c(1, 2, 3, 6)) %>%
            pivot_wider(names_from = "s", values_from = c("y_hat")) %>% 
            rename(y_hat_1 = `1`, y_hat_2 = `2`, y_hat_3 = `3`, y_hat_6 = `6`) %>%
            mutate( y_hat_3 = y_hat_3 + y_hat_6, .keep = "unused") %>% 
            pivot_longer(y_hat_1:y_hat_3, names_to = "state", values_to = "value") %>%
            mutate(state = recode(state, !!!state_names)) %>% mutate(antigen = "Spike model")

    y_hat_dist_rbd <- stanfit_rbd$draws(c("y_hat")) %>% spread_draws(y_hat[n, t, s]) %>% rename(PID = n, t_id = t) %>%
        mutate(PID = factor(PID, levels = 1:41)) %>%  filter(s %in% c(1, 2, 3, 6)) %>%
            pivot_wider(names_from = "s", values_from = c("y_hat")) %>% 
            rename(y_hat_1 = `1`, y_hat_2 = `2`, y_hat_3 = `3`, y_hat_6 = `6`) %>%
            mutate( y_hat_3 = y_hat_3 + y_hat_6, .keep = "unused") %>% 
            pivot_longer(y_hat_1:y_hat_3, names_to = "state", values_to = "value") %>%
            mutate(state = recode(state, !!!state_names)) %>% mutate(antigen = "RBD model")
    y_hat_dist <- bind_rows(y_hat_dist_s, y_hat_dist_rbd) 
    data_i_all %>% head

    y_hat_dist_data <- data_i_all %>% left_join(y_hat_dist %>% filter(!is.na(value))) %>% ungroup %>% 
        mutate(key = paste(PID, t_id, sep = "_")) %>%
        mutate(key = factor(key, levels = unique(key))) %>%
        mutate(key = as.numeric(key)) 

    df_rmse_all <- map_df(as.character(state_names),
        function(state_str) {
            #  "Memory B-cell conc."               "Plasmablast conc." "sVNT to ancestral variant (log)" 
           # bcell_mat <- y_hat_dist_data %>% ungroup %>% filter(!is.na(i), !is.na(data_value)) %>% select(.draw, PID, i, antigen, day_data, state, value, data_value) %>%
           #  state_str <- "sVNT to ancestral variant (log)"
            y_hat_dist_data_state <- y_hat_dist_data %>% filter(state == state_str) 
                
            mat_model_post <- y_hat_dist_data_state %>% select(antigen, key, PID, t, value) %>% pull(value) %>% matrix(ncol = 4000, byrow = TRUE) %>% t
            df_model_data <- y_hat_dist_data_state %>% select(antigen, key, PID, t, data_value) %>% unique 
            vec_model_data <- df_model_data %>% pull(data_value)

            pred_mean <- colMeans(mat_model_post)
            rmse_score <- sqrt((vec_model_data - pred_mean)^2)

            df_model_data$rmse_score <- rmse_score
            df_model_post_summary <- y_hat_dist_data_state %>% select(key, PID, t, value) %>% group_by(PID, key, t) %>% mean_qi(value)

            left_join(df_model_post_summary, df_model_data) %>% mutate(state = state_str)
        }
    )  %>% mutate(VisitN = cut(t, breaks = c(-1, 13, 22, 300), labels = c("2–13 days", "14–21 days", "67–210 days") ) ) %>%
         mutate(VisitN = factor(VisitN, levels = c("2–13 days", "14–21 days", "67–210 days")) ) %>% filter(!is.na(VisitN)) %>%
        mutate(antigen = factor(antigen, levels = c("Spike model", "RBD model"))) 

    df_rmse_all %>% 
        ggplot() + 
            geom_boxplot(aes(x = state, y = rmse_score, fill = antigen), position = position_dodge(0.8), alpha = 0.3) + 
            theme_bw() +
            scale_fill_manual(values = c(color_study_A, color_study_B)) + 
            david_theme_B() + 
            theme(legend.position = "top") + 
            labs(x = "", y = "Root Mean Squared Error\n (RMSE)", fill = "", color = "") +        
             theme(strip.text.x = element_text(size = 24), axis.text = element_text(size = 20), axis.title.y = element_text(size = 22), 
              title = element_text(size = 20), legend.text = element_text(size = 20))   
    ggsave(here::here("outputs", "figs", "post_pred", "fig1B_RMSE.png"), height = 10, width = 15 )


    p1 <- df_rmse_all %>% filter(antigen == "Spike model") %>%
            ggplot() + 
                geom_linerange(aes(PID, y = value, ymin = .lower, ymax = .upper), size = 1) + 
                geom_point(aes(PID, y = value), color = "black", size = 3)  +
                geom_point(aes(PID, y = data_value, fill = rmse_score), alpha = 0.7, shape = 22, color = "white", size = 3)  +
                facet_nested(VisitN ~ state, scales = "free") + theme_bw() +
                labs(x = "Person ID", y = "Value", fill = "Root Mean Squared Error (RMSE)", color = "") + 
                scale_color_manual(labels = c( "Model fit"), values = c("black") ) + 
                # have low values red and high values gray in gradient
                scale_fill_continuous( low = "red", high = "gray")+ 
                theme(legend.position = "bottom") + 
                ggtitle(paste0("Model fits to calibration data for " , "Ancestral spike",  " antigen") ) + 
                david_theme_A() + theme(axis.text.x = element_blank())
    p2 <- df_rmse_all %>% filter(antigen == "RBD model") %>%
            ggplot() + 
                geom_linerange(aes(PID, y = value, ymin = .lower, ymax = .upper), size = 1) + 
                geom_point(aes(PID, y = value), color = "black", size = 3)  +
                geom_point(aes(PID, y = data_value, fill = rmse_score), alpha = 0.7, shape = 22, color = "white", size = 3)  +
                facet_nested(VisitN ~ state, scales = "free") + theme_bw() +
                labs(x = "Person ID", y = "Value", fill = "Root Mean Squared Error (RMSE)", color = "") + 
                scale_color_manual(labels = c( "Model fit"), values = c("black") ) + 
                # have low values red and high values gray in gradient
                scale_fill_continuous( low = "red", high = "gray")+ 
                theme(legend.position = "bottom") + 
                ggtitle(paste0("Model fits for calibration data" , " Ancestral RBD",  " antigen") ) + 
                david_theme_A() + theme(axis.text.x = element_blank())
    p1 / p2
    ggsave(here::here("outputs", "figs", "post_pred", "figS1_RMSE.png"), height = 15, width = 15)

}


#' Plot the individual-level trajectories
plot_cal_traj_ind <- function() {

   # Get the trajectories for each individual
    df_uncert_ind_s <- readRDS(here::here("outputs", "stanfit",  paste0("traj_ind_",  "s", ".RDS"))) %>% mutate(antigen = "Ancestral spike")
    df_uncert_ind_rbd <- readRDS(here::here("outputs", "stanfit",  paste0("traj_ind_",  "rbd", ".RDS"))) %>% mutate(antigen = "Ancestral RBD")
    df_uncert_ind <- bind_rows(df_uncert_ind_s, df_uncert_ind_rbd)  %>%
        mutate(antigen = factor(antigen, levels = c("Ancestral spike", "Ancestral RBD"))) 


    # Get the calibration data
    data_i_all_s <- get_data_ind("nih_wu_s") %>% mutate(antigen = "Ancestral spike")
    data_i_all_rbd <- get_data_ind("nih_wu_rbd") %>% mutate(antigen = "Ancestral RBD")
    data_i_all <- bind_rows(data_i_all_s, data_i_all_rbd)  %>%
        mutate(antigen = factor(antigen, levels = c("Ancestral spike", "Ancestral RBD"))) 


    p1 <- df_uncert_ind %>% filter(antigen == "Ancestral spike", time < 220) %>% 
        mutate(PID = factor(PID, levels = unique(PID))) %>%
        mutate(PID = as.numeric(PID)) %>% 
        mutate(PID = factor(PID, levels = as.character(1:41))) %>%
        ggplot() +
        geom_ribbon(aes(x = time, y = value, ymin = .lower, ymax = .upper, fill = state_variable), 
        alpha = 0.3) + 
        geom_line(aes(x = time, y = value, color = state_variable)) + 
        geom_point(data = data_i_all  %>% filter(antigen == "Ancestral spike") %>% mutate(PID = factor(PID, levels = unique(PID))) %>%
        mutate(PID = as.numeric(PID)) %>% mutate(PID = factor(PID, levels = as.character(1:41))), aes(x = t, y = value, fill = state), shape = 21, size = 3, alpha = 0.5) + theme_bw() +    
        xlim(0, 220) +
        facet_wrap(vars(PID)) + labs(x = "Days post-vaccination", y = "Conc.(%)/sVNT (log)", fill = "Data", color = "Posterior predictive interval") + 
        theme(legend.position = "top") + 
        ggtitle(paste0("Individual-level immune trajectories in response to ",  "Ancestral spike", " antigen on calibration dataset") )

    p2 <- df_uncert_ind %>% filter(antigen == "Ancestral RBD", time < 220) %>% mutate(PID = factor(PID, levels = unique(PID))) %>%
        mutate(PID = as.numeric(PID)) %>% mutate(PID = factor(PID, levels = as.character(1:41))) %>%
        ggplot() +
        geom_ribbon(aes(x = time, y = value, ymin = .lower, ymax = .upper, fill = state_variable), 
        alpha = 0.3) + 
        geom_line(aes(x = time, y = value, color = state_variable)) + 
        geom_point(data = data_i_all %>% filter(antigen == "Ancestral RBD") %>% mutate(PID = factor(PID, levels = unique(PID))) %>%
        mutate(PID = as.numeric(PID)) %>% mutate(PID = factor(PID, levels = as.character(1:41))), aes(x = t, y = value, fill = state), shape = 21, size = 3, alpha = 0.5) + theme_bw() +    
        xlim(0, 220) +
        facet_wrap(vars(PID)) + labs(x = "Days post-vaccination", y = "Conc.(%)/sVNT (log)", fill = "Data", color = "Posterior predictive interval") + 
        theme(legend.position = "top") + 
        ggtitle(paste0("Individual-level immune trajectories in response to ",  "Ancestral RBD", " antigen on calibration dataset") )
   # ggsave(here::here("outputs", "figs", "manu_A", foldername, "fits", "fig1_ind_traj_XX.png"), width = 15, height = 12, units = "in")
    #p1 <- p4A  + plot_layout(heights = c(1, 2)) + plot_annotation(tag_levels = "A") + plot_layout(guides = "collect") &  theme(legend.position = "top") 
    p1 / p2 + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") & theme(legend.position='top') 
    ggsave(here::here("outputs", "figs", "post_pred", "figS3.png"), height = 20)
}