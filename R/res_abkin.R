
# This script contains the functions to plot the antibody kinetics and the transition time for the antibody production following vaccination with the second dose.
# The functions are used in the main script `main.R` to generate the figures for the manuscript.
#' Code to plot Figure 4 of manuscript, the antibody kinetics

#' @param stanfit The stanfit object
#' @param fig_folder The folder where the data is stored
#' @export
plot_ab_kinetics <- function(stanfit, fig_folder) {

    if (fig_folder == "nih_wu_s") {
        antigen_str <- "Ancestral spike"
        antigen_file <- "s"
    } else {
        antigen_str <- "Ancestral RBD"
        antigen_file <- "rbd"
    }

    relabel_data_vaccine <- c("P" = "BNT162b2", "AZ" = "ChAdOx1")

    relabel_data_state <- c("Bmem" = "Memory B-cell conc.",
        "pbBmem" = "Plasmablast conc.", "IC50hs" = "sVNT to ancestral variant")
    
    data_model <- readRDS(here::here("outputs", "data_clean", fig_folder, "df_meta.RDS"))
    data_fit <- readRDS(here::here("outputs", "data_clean", fig_folder, "df_model_data.RDS"))
    data_fit_clean <- data_fit %>% select(PID:age_cat) %>%
        pivot_longer(c(Bmem, pbBmem, IC50hs), names_to = "state", values_to = "value")  %>%
        mutate(state = recode(state, !!!relabel_data_state)) %>% filter(value > -1) %>%
        mutate(vax_type = recode(vax_type, !!!relabel_data_vaccine)) 
    
    relabel_model_vaccine <- c("2" = "ChAdOx1", "1" = "BNT162b2")
    relabel_l_s <-c("1" = "Memory B-cell conc.",
        "2" = "Plasmablast conc.", 
        "3" = "sVNT from plasmablasts", "4" = "GC", "5" = "PB", 
        "6" = "sVNT from plasma cells", "tot" = "sVNT total")

    full_dynamics_vac_2_wu_s <- stanfit$draws(c("y_full_vac")) %>%
            spread_draws(y_full_vac[i, t, s])
    full_dynamics_vac_2a_wu_s <- full_dynamics_vac_2_wu_s %>% filter(s %in% c(1, 2, 3, 6)) %>%
        pivot_wider(names_from = "s", values_from = "y_full_vac") %>%
        mutate(tot = `3` + `6`, prop = `3` / (`3` + `6`)) %>% 
        pivot_longer(`1`:`prop`, names_to = "s", values_to = "y_full_vac") %>%
        group_by(i, t, s) %>%
        mean_qi(y_full_vac) 

    tt_df <- full_dynamics_vac_2a_wu_s %>% mutate(s = recode(s, !!!relabel_l_s)) %>%
        filter(s == "prop") %>%
        mutate(vax_type = recode(i, !!!relabel_model_vaccine)) %>% drop_na %>% filter(y_full_vac < 0.5) %>% 
        group_by(vax_type) %>% filter(t == min(t)) 

    p1 <- full_dynamics_vac_2a_wu_s %>% mutate(s = recode(s, !!!relabel_l_s)) %>%
        filter(s != "prop") %>%
        filter(!s %in% c("Titre (short-lived) (log)", "Titre (long-lived) (log)")) %>% 
        mutate(vax_type = recode(i, !!!relabel_model_vaccine)) %>% drop_na %>% 
        filter(s %in% c(
            "sVNT from plasmablasts",
            "sVNT from plasma cells", 
            "sVNT total")) %>%
        ggplot() + 
        geom_vline(data = tt_df, aes(xintercept = t), linetype = "dashed", color = "gray20", size = 1) +
        geom_text(data = tt_df, aes(x = t + 10, y = 2.3), label = "Transition time (TT)", angle = 90, color = "gray20", size = 6) +
        geom_ribbon(aes(x = t, ymin = .lower, 
            ymax = .upper, fill = s), alpha = 0.2) + 
        geom_line(aes(x = t, y = y_full_vac, group = s, color = s, alpha = s), 
            size = 2) + 
        guides(fill = "none", alpha = "none") +
        scale_color_manual(values = c("#f6afad", "#960319", "#552531", "#552531")) + 
        scale_fill_manual(values = c("#f6afad", "#960319", "#552531", "#552531")) + 
        scale_alpha_manual(values = c(0.5, 0.5, 1)) + 
        geom_point(data = 
            data_fit_clean %>% filter(state == "sVNT to ancestral variant"), 
            aes(time_post_exp, value), 
            alpha = 0.75, size = 4, shape = 21, fill = "#552531") +
        facet_grid(cols = vars(vax_type)) + theme_bw() + 
        theme(strip.text.x = element_text(size = 24), axis.text = element_text(size = 15), axis.title = element_text(size = 18), 
              title = element_text(size = 20), legend.text = element_text(size = 15))  +
        labs(x = "Time after vaccination with second dose (days)", y = "sVNT (log)", color = "Model-trajectories", fill = NA) + 
        ggtitle("Temporal variation in origin of antibody production following vaccination with second dose")


    df_uncert_all_wu_t <- clean_fig4_alt_time(stanfit, fig_folder)
    df_uncert_all_wu_a <- clean_fig4_alt_age(stanfit, fig_folder)

    ### HERE LOOKS GOOD ###

    pg_time_s <- df_uncert_all_wu_t %>% convert_to_data_time
    pg_age_s <- df_uncert_all_wu_a %>% convert_to_data_age

    recode_time <- c("<28 days" = "<28 days", ">= 28 days" = "28+ days")

    p2 <- pg_time_s %>% 
        mutate(time_d1 = recode(time_d1, !!!recode_time)) %>%
        ggplot(aes(x = x, y = y, color = as.character(time_d1) )) + 
        geom_hline(yintercept = c(0.025, 0.975), linetype = "dotted") + 
        geom_hline(yintercept = 0.5, linetype = "dashed") + 
        geom_step(size = 2, alpha = 0.7) + 
        labs(x = "Time after vaccination with second dose (days)", 
            y = "Proportion of posterior \ndistribution above threshold", color = "Time since first dose") +
        facet_grid(cols = vars(v)) + theme_bw() + ylim(0, 1) + xlim(0, 365) +
        scale_color_manual(values = c("#ced011", "#b942bb")) + 
        theme(strip.text.x = element_text(size = 24), axis.text = element_text(size = 15), axis.title = element_text(size = 18), 
              title = element_text(size = 20), legend.text = element_text(size = 15)) + 
        labs(title = "Persistence of antibodies levels following vaccination with second dose", subtitle =  "Threshold: IC50 value of 10") 
              


    p3 <- pg_age_s %>% 
        mutate(age_group = paste0(age_group, " yrs")) %>%
        ggplot(aes(x = x, y = y, color = age_group )) + 
        geom_hline(yintercept = c(0.025, 0.975), linetype = "dotted") + 
        geom_hline(yintercept = 0.5, linetype = "dashed") + 
        geom_step(size = 2, alpha = 0.7) + 
        labs(x = "Time after vaccination with second dose (days)", 
            y = "Proportion of posterior \ndistribution above threshold", color = "Age group") +
        facet_grid(cols = vars(v)) + theme_bw() + ylim(0, 1) + xlim(0, 365) +
        scale_color_manual(values = c("#960319","#c0712e", "#ced011",  "#89c462", "#0b2b53")) + 
        theme(strip.text.x = element_text(size = 24), axis.text = element_text(size = 15), axis.title = element_text(size = 18), 
        title = element_text(size = 22), legend.text = element_text(size = 15)) 

    p1 / p2 / p3 + plot_annotation(tag_levels = "A", title = antigen_str, theme = theme(title = element_text(size = 20))) + plot_layout(heights = c(2, 1, 1))
    ggsave(here::here("outputs", "figs", "abkin", paste0("ab_compare_final", antigen_file, ".pdf")),  height = 16, width = 16)


}


#' A function to simulate the in-host kinetic model outside of stan
#' @param stanfit The stanfit object
#' @param fig_folder The folder where the data is stored
#' 
clean_fig4_alt <- function(stanfit, fig_folder) {

   # stanfit_name <- "nih_vac_wu_s_i"

    #post_sample_cross <- readRDS(file = here::here("outputs", "stanfit", "manu_A", paste0(stanfit_name, ".RDS")) )
    data_model <- readRDS(here::here("outputs", "data_clean", fig_folder, "df_meta.RDS"))
    data_listA <- readRDS(here::here("outputs", "data_clean", fig_folder, "list_stan_data.RDS"))

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
    df_uncert_all <- map_df(1:2, 
        function(s) {
            init_i <- data_listA$y0_init_vac[(s - 1) %% 2 + 1, ]
            names(init_i) <- c("B", "P", "A1", "G", "L", "A2")

            param_uncert_sum_i <- param_uncert %>% filter( v == param_uncert_sum$v[s])

            df_uncert <- map_df(1:M,
                ~ode_results_vac(param_uncert_sum_i[.x, 5:15], 365, init_i) %>% as.data.frame  )  %>%
                select(time, B, P, A1, A2) %>% mutate(A = A1 + A2, .keep = "unused") %>% 
                mutate(i = map(1:M, ~rep(.x, 366)) %>% unlist) %>% 
            filter(A > log10(10)) %>% group_by(i) %>% filter(time == max(time)) %>% ungroup %>% complete(i = 1:M, fill = list(time = 0)) %>% 
            mutate( v = param_uncert_sum$v[s])
        }
    )
    df_uncert_all
}

#' Function to do the data manipulation and simulate the trajectories for figure 4 stratified by time
#' @param stanfit The stanfit object
#' @param fig_folder The folder where the data is stored
#' 
clean_fig4_alt_time <- function(stanfit, fig_folder) {

   # stanfit_name <- "nih_vac_wu_s_i"
   # fig_folder <- "nih_wu_s"

    #post_sample_cross <- readRDS(file = here::here("outputs", "stanfit", "manu_A", paste0(stanfit_name, ".RDS")) )
    data_model <- readRDS(here::here("outputs", "data_clean", fig_folder, "df_meta.RDS"))
    data_listA <- readRDS(here::here("outputs", "data_clean", fig_folder, "list_stan_data.RDS"))

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
    df_uncert_all <- map_df(1:4, 
        function(s) {
            init_i <- data_listA$y0_init_time[(s - 1) %% 2 + 1, ]
            names(init_i) <- c("B", "P", "A1", "G", "L", "A2")

            param_uncert_sum_i <- param_uncert %>% filter(t == param_uncert_sum$t[s], v == param_uncert_sum$v[s])

            df_uncert <- map_df(1:M,
                ~ode_results_vac(param_uncert_sum_i[.x, 6:16], 365, init_i) %>% as.data.frame  )  %>%
                select(time, B, P, A1, A2) %>% mutate(A = A1 + A2, .keep = "unused") %>% 
                mutate(i = map(1:M, ~rep(.x, 366)) %>% unlist) %>% 
            filter(A > log10(10)) %>% group_by(i) %>% filter(time == max(time)) %>% ungroup %>% complete(i = 1:M, fill = list(time = 0)) %>% 
            mutate(t = param_uncert_sum$t[s], v = param_uncert_sum$v[s])
        }
    )
    df_uncert_all
}

#' Function to do the data manipulation and simulate the trajectories figure 4 stratified by age
#' @param stanfit The stanfit object
#' @param fig_folder The folder where the data is stored
#' 
clean_fig4_alt_age <- function(stanfit, fig_folder) {
    data_model <- readRDS(here::here("outputs", "data_clean", fig_folder, "df_meta.RDS"))
    data_listA <- readRDS(here::here("outputs", "data_clean", fig_folder, "list_stan_data.RDS"))

    param_uncert <- stanfit$draws(
        c("theta_a",   
            "prop_asc", "deltaAg", "deltaP", "deltaL", "deltaA", "deltaG")) %>%
        spread_draws(theta_a[j, i, s], prop_asc, deltaAg, deltaP, deltaL, deltaA, deltaG) %>% 
        mutate(deltaM = 0.001, .before = "deltaP") %>%
        mutate(deltaAg = 1 / (deltaAg * 30), deltaP = 1 / deltaP, deltaL = 1 / deltaL, deltaA = 1 / deltaA, deltaG = 1 / deltaG) %>% 
        pivot_wider(names_from = "s", values_from = "theta_a") %>% 
        select(j:.draw, `1`:`4`, prop_asc, deltaAg:deltaG) %>% 
        rename(a = j, v = i) %>% rename(param1 = `1`, param2 = `2`, param3 = `3`, param4 = `4`)

    param_uncert_sum <- param_uncert %>% ungroup %>% summarise(across(param1:deltaG, ~mean(.x)), .by = c("a", "v")) 

    M <- 100
    df_uncert_all <- map_df(1:10, 
        function(s) {
            init_i <- data_listA$y0_init_age[(s + 1) %/% 2, ]
            names(init_i) <- c("B", "P", "A1", "G", "L", "A2")

            param_uncert_sum_i <- param_uncert %>% filter(a == param_uncert_sum$a[s], v == param_uncert_sum$v[s])

            df_uncert <- map_df(1:M,
                ~ode_results_vac(param_uncert_sum_i[.x, 6:16], 365, init_i) %>% as.data.frame  )  %>%
                select(time, B, P, A1, A2) %>% mutate(A = A1 + A2, .keep = "unused") %>% 
                mutate(i = map(1:M, ~rep(.x, 366)) %>% unlist) %>% 
            filter(A > log10(10)) %>% group_by(i) %>% filter(time == max(time)) %>% ungroup %>% complete(i = 1:M, fill = list(time = 0)) %>% 
            mutate(a = param_uncert_sum$a[s], v = param_uncert_sum$v[s])
        }
    )
    df_uncert_all
}

#' Code to plot Figure 4 of manuscript, the antibody kinetics
#' 
plot_abkin_tt <- function()  {

    if (fig_folder == "nih_wu_s") {
        antigen_str <- "Ancestral spike"
        antigen_file <- "s"
    } else {
        antigen_str <- "Ancestral RBD"
        antigen_file <- "rbd"
    }

    relabel_data_vaccine <- c("P" = "BNT162b2", "AZ" = "ChAdOx1")

    relabel_data_state <- c("Bmem" = "Memory B-cell conc.",
        "pbBmem" = "Plasmablast conc.", "IC50hs" = "sVNT to ancestral variant")
    
    data_model <- readRDS(here::here("outputs", "data_clean", fig_folder, "df_meta.RDS"))
    data_fit <- readRDS(here::here("outputs", "data_clean", fig_folder, "df_model_data.RDS"))
    data_fit_clean <- data_fit %>% select(PID:age_cat) %>%
        pivot_longer(c(Bmem, pbBmem, IC50hs), names_to = "state", values_to = "value")  %>%
        mutate(state = recode(state, !!!relabel_data_state)) %>% filter(value > -1) %>%
        mutate(vax_type = recode(vax_type, !!!relabel_data_vaccine)) 
    
    relabel_model_vaccine <- c("2" = "ChAdOx1", "1" = "BNT162b2")
    relabel_l_s <-c("1" = "Memory B-cell conc.",
        "2" = "Plasmablast conc.", 
        "3" = "sVNT from plasmablasts", "4" = "GC", "5" = "PB", 
        "6" = "sVNT from plasma cells", "tot" = "sVNT total")

    stanfit_s <- readRDS(file = here::here("outputs", "stanfit", paste0( "nih_vac_wu_s_i", ".RDS")) )
    full_dynamics_vac_2_wu_s <- stanfit_s$draws(c("y_full_vac")) %>%
            spread_draws(y_full_vac[i, t, s])
    full_dynamics_vac_2a_wu_s <- full_dynamics_vac_2_wu_s %>% filter(s %in% c(1, 2, 3, 6)) %>%
        pivot_wider(names_from = "s", values_from = "y_full_vac") %>%
        mutate(tot = `3` + `6`, prop = `3` / (`3` + `6`)) %>% 
        pivot_longer(`1`:`prop`, names_to = "s", values_to = "y_full_vac") %>%
        group_by(i, t, s) %>%
        mean_qi(y_full_vac) %>% mutate(antigen = "Ancestral spike")

    stanfit_rbd <- readRDS(file = here::here("outputs", "stanfit", paste0( "nih_vac_wu_rbd_i", ".RDS")) )
    full_dynamics_vac_2_wu_rbd <- stanfit_rbd$draws(c("y_full_vac")) %>%
            spread_draws(y_full_vac[i, t, s])
    full_dynamics_vac_2a_wu_rbd <- full_dynamics_vac_2_wu_rbd %>% filter(s %in% c(1, 2, 3, 6)) %>%
        pivot_wider(names_from = "s", values_from = "y_full_vac") %>%
        mutate(tot = `3` + `6`, prop = `3` / (`3` + `6`)) %>% 
        pivot_longer(`1`:`prop`, names_to = "s", values_to = "y_full_vac") %>%
        group_by(i, t, s) %>%
        mean_qi(y_full_vac)  %>% mutate(antigen = "Ancestral RBD")

    full_dynamics_vac_2a_wu <- full_dynamics_vac_2a_wu_s %>% bind_rows(full_dynamics_vac_2a_wu_rbd) %>% 
        mutate(antigen = factor(antigen, levels = c("Ancestral spike", "Ancestral RBD") ))

    
    p2 <- full_dynamics_vac_2a_wu %>% mutate(s = recode(s, !!!relabel_l_s)) %>%
        filter(s == "prop") %>%
        mutate(vax_type = recode(i, !!!relabel_model_vaccine)) %>% drop_na %>% ggplot() + 
        geom_ribbon(aes(x = t, ymin = .lower, 
            ymax = .upper, fill = vax_type), alpha = 0.3) + 
        geom_line(aes(x = t, y = y_full_vac, color = vax_type), size = 2) + 
        theme_bw() + 
        guides(fill = "none") + 
        scale_color_manual(values = c("#2271B2", "#D55E00") ) +
        labs(x = "Days post vaccination",
            y = "Proportion of antibodies secreted from plasmablasts",
        color = "Vaccine type") +
        theme(strip.text.x = element_text(size = 24), axis.text = element_text(size = 15), axis.title = element_text(size = 18), 
              title = element_text(size = 20), legend.text = element_text(size = 15)) + 
        facet_grid(cols = vars(antigen)) 

    #"Temporal variation in origin of antibody production following vaccination with second dose"
    # Persistence of antibodies levels following vaccination with second dose
    (p2) + plot_annotation(tag_levels = "A") + plot_layout( )
    ggsave(here::here("outputs", "figs", "abkin", paste0("ab_compare_tt.pdf")),  height = 16, width = 16)

}


# Helper functions
#' Function to convert the CDF to manipulatable data, stratified by time
#' @param data_l The data list
convert_to_data_time <- function(data_l) {
    p  <- ggplot(data_l, aes(time, color = as.character(t))) + stat_ecdf() + facet_grid(vars(v))
    pg <- ggplot_build(p)$data[[1]] %>% rename(time_d1 = group, v = PANEL) %>% filter(x >= 0, x < 365) %>%
        select(y, x, v, time_d1) %>% mutate(y = 1-y) %>%  mutate(v = recode(v, "2" = "ChAdOx1", "1" = "BNT162b2")) %>% mutate(time_d1 = recode(time_d1, "1" = "<28 days", "2" = ">= 28 days"))
    to_add_max <- pg %>% group_by(v, time_d1) %>% summarise(y = max(y)) %>% mutate(x = 0)
    to_add_min <- pg %>% group_by(v, time_d1) %>% summarise(y = min(y)) %>% mutate(x = 365)
    pg <- bind_rows(pg, to_add_min) %>% bind_rows(to_add_max)
}

#' Helper functions
#' Function to convert the CDF to manipulatable data, stratified by age
#' @param data_l The data list
#' 
convert_to_data_age <- function(data_l) {
    p  <- ggplot(data_l, aes(time, color = as.character(a))) + stat_ecdf() + facet_grid(vars(v))
    pg <- ggplot_build(p)$data[[1]] %>% rename(age_group = group, v = PANEL) %>% filter(x >= 0, x < 365) %>%
        select(y, x, v, age_group) %>% mutate(y = 1-y) %>%  mutate(v = recode(v, "2" = "ChAdOx1", "1" = "BNT162b2")) %>% mutate(age_group = recode(age_group, "1" = "<30", "2" = "30–39", "3" = "40–49", "4" = "50-59", "5" = "60+"))
    to_add_max <- pg %>% group_by(v, age_group) %>% summarise(y = max(y)) %>% mutate(x = 0)
    to_add_min <- pg %>% group_by(v, age_group) %>% summarise(y = min(y)) %>% mutate(x = 365)
    pg <- bind_rows(pg, to_add_min) %>% bind_rows(to_add_max)
}
