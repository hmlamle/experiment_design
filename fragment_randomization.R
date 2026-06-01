# Hannah-Marie Lamle
# Coral Fragment Randomization for Experiment
# May 31, 2026


library(tidyverse)
library(readxl)
library(ompr)
library(ompr.roi)
library(ROI.plugin.glpk)


# -------------------- Load Datasheet for frags ------------------------

# load master dataset and clean it so that it's just sids and only those still alive:
sids_master <- read_excel("ssid_collection_MASTER.xlsx")

# clean up the date column from weird excel format:
sids_master$date[sids_master$date == "?"] <- NA
sids_master$date <- sids_master$date %>%
  as.character() %>%
  trimws() %>%                     # removes hidden spaces
  gsub("[^0-9]", "", .) %>%       # keeps ONLY numbers
  as.numeric()

sids_master$date <- as.Date(sids_master$date, origin = "1899-12-30")

# remove frags i know i dont want in the experiment:
ssids <- sids_master %>%
  filter(species == "ssid") %>%                                             # removing the radians
  filter(is.na(`intake_reflectance?`) | `intake_reflectance?` !="F") %>%    # removing frags already frozen
  subset(fragment_ID != "sid-dead")                                         # remove random sid I got
  

## Further cleaning of the dataframe: 

# consolidating all the POMs and PEVs

ssids <- ssids %>%
  mutate(site = case_match(site,
                           "POM 7"  ~ "POM",
                           "POM 1"  ~ "POM",
                           "POM3"   ~ "POM",
                           "PEV 13" ~ "PEV",
                           "PEV 5"  ~ "PEV",
                           "FTL-B1" ~ "BAR",
                           "FTL-B2" ~ "BAR",
                           "MB 2"   ~ "MB2",
                           .default = site)) %>%
  filter(site != "MB2")  # and removing the miami beach site since only 2 were sids

# After talking with Ryan, we should remove also the corals 
# that we collected a long time ago: 

ssids <- ssids %>%
  filter(date >= as.Date("2026-03-22")) %>%
  filter(!site %in% c("POM", "PEV"))

# Create DF's for tank treatments:

tank <- c(1,2,3,4,5,6,7,8,9,10,11,12)
treatments <- c("+4C", "control", "+2C", "+4C", "control", "+2C", "control", "+4C", "+2C", "control", "+4C", "+2C")

tank_treatments <- data.frame(tank, treatments)


# ----------------------- Balanced Allocation algorithm ----------------

set.seed(123)

n_tanks <- 12
tank_size <- 25
target_total <- n_tanks * tank_size
tank_ids <- paste0("Tank_", 1:n_tanks)

## ---------- 1. Force Dodge Island: 3 per tank --------------------

dodge <- ssids %>%
  filter(site == "Dodge Island") %>%
  slice_sample(n = 36) %>%
  mutate(
    tank = rep(tank_ids, each = 3)
  )
## ---------- 2. Remove Dodge Island -----------------------------

ssids_no_dodge <- ssids %>%
  filter(site != "Dodge Island")

# Need 300 total, already using 36 Dodge Island
remaining_total <- target_total - nrow(dodge)

## --------- 3. Decide how many fragments to keep per remaining site ------------

site_counts <- ssids_no_dodge %>%
  count(site, name = "available") %>%
  mutate(keep_n = available)

to_drop <- sum(site_counts$keep_n) - remaining_total
site_counts <- site_counts %>%
  arrange(desc(available)) %>%
  mutate(
    max_drop = pmax(available - n_tanks, 0),
    drop_n = pmin(max_drop, pmax(to_drop - cumsum(lag(max_drop, default = 0)), 0)),
    keep_n = available - drop_n
  ) %>%
  select(site, keep_n)

## -------------- 4. Subsample remaining sites -----------------------------

ssids_keep <- ssids_no_dodge %>%
  left_join(site_counts, by = "site") %>%
  group_by(site) %>%
  slice_sample(prop = 1) %>%
  mutate(site_row = row_number()) %>%
  filter(site_row <= keep_n) %>%
  ungroup() %>%
  select(-site_row, -keep_n)

## ----------- 5. Assign remaining fragments evenly across tanks -----------

tank_load <- rep(3, n_tanks)   # Dodge Island already has 3 per tank
names(tank_load) <- tank_ids
assigned_list <- ssids_keep %>%
  group_split(site)
out <- list()

for(i in seq_along(assigned_list)) {
  dat <- assigned_list[[i]]
  n_site <- nrow(dat)
  base <- floor(n_site / n_tanks)
  extra <- n_site %% n_tanks
  site_tank_counts <- rep(base, n_tanks)
  if(extra > 0) {
    extra_tanks <- order(tank_load)[1:extra]
    site_tank_counts[extra_tanks] <- site_tank_counts[extra_tanks] + 1
  }
  tank_vec <- rep(tank_ids, site_tank_counts)
  tank_vec <- sample(tank_vec)
  dat$tank <- tank_vec
  tank_load <- tank_load + site_tank_counts
  names(tank_load) <- tank_ids
  out[[i]] <- dat
}

other_assigned <- bind_rows(out)

## ------------ 6. Combine Dodge Island + all others ------------------------

ssids_assigned <- bind_rows(dodge, other_assigned)

## -------------- 7. Checks -----------------------------

# Should be exactly 25 per tank
ssids_assigned %>%
  count(tank)

# Dodge Island should be exactly 3 per tank
ssids_assigned %>%
  filter(site == "Dodge Island") %>%
  count(tank)

# Total fragments used
nrow(ssids_assigned)


ssids_assigned %>%
  
  count(site)

table(ssids_assigned$tank, ssids_assigned$site)


### change tank column to work with my tank assignment df: 

ssids_assigned$tank <- as.integer(gsub("Tank_", "", ssids_assigned$tank))



# ------------------ Assigning time point sampling ------------------
# add the treatment for each tank into the DF: 
ssids_assigned <- ssids_assigned %>%
  left_join(tank_treatments, by = "tank") %>%
  mutate(treatment = treatments) %>%
  select(-treatments)

# change the name of the sample timepoint column so it's easier in R:
ssids_assigned <- ssids_assigned %>%
  rename(sample_time = "sample timepoint")





set.seed(123)
n_timepoints <- 5
timepoints <- paste0("T", 1:n_timepoints)
ssids_sampling <- ssids_assigned %>%
  mutate(
    tank = as.numeric(tank),
    sample_time = NA_character_,
    row_id = row_number()
  )
out <- list()
for(this_tank in sort(unique(ssids_sampling$tank))) {
  dat_tank <- ssids_sampling %>%
    filter(tank == this_tank)
  tank_load <- rep(0, n_timepoints)
  names(tank_load) <- timepoints
  site_order <- dat_tank %>%
    count(site, name = "n") %>%
    arrange(desc(n))
  for(s in site_order$site) {
    rows <- dat_tank %>%
      filter(site == s) %>%
      slice_sample(prop = 1) %>%
      pull(row_id)
    n_site <- length(rows)
    base <- floor(n_site / n_timepoints)
    extra <- n_site %% n_timepoints
    target_counts <- rep(base, n_timepoints)
    if(extra > 0) {
      target_counts[order(tank_load)[1:extra]] <- 
        target_counts[order(tank_load)[1:extra]] + 1
    }
    names(target_counts) <- timepoints
    current_counts <- rep(0, n_timepoints)
    names(current_counts) <- timepoints
    for(r in rows) {
      possible_times <- timepoints[
        tank_load < 5 &
          current_counts < target_counts
      ]
      if(length(possible_times) == 0) {
        possible_times <- timepoints[tank_load < 5]
      }
      best_time <- possible_times[
        which.min(tank_load[possible_times])
      ]
      dat_tank$sample_time[dat_tank$row_id == r] <- best_time
      tank_load[best_time] <- tank_load[best_time] + 1
      current_counts[best_time] <- current_counts[best_time] + 1
    }
  }
  out[[as.character(this_tank)]] <- dat_tank
}
ssids_sampling <- bind_rows(out) %>%
  select(-row_id)

ssids_assigned = ssids_sampling




set.seed(123)

timepoints <- paste0("T", 1:5)
need_fix <- ssids_assigned %>%
  count(site, treatment, name = "total_n") %>%
  filter(total_n > 0)

for(i in seq_len(nrow(need_fix))) {
  this_site <- need_fix$site[i]
  this_treatment <- need_fix$treatment[i]
  repeat {
    group_rows <- which(
      ssids_assigned$site == this_site &
        ssids_assigned$treatment == this_treatment
    )
    tab <- table(factor(
      ssids_assigned$sample_time[group_rows],
      levels = timepoints
    ))
    # stop if it is already 1,1,1,1,1
    if(all(tab == 1)) break
    over_time <- names(tab)[tab > 1][1]
    under_time <- names(tab)[tab == 0][1]
    # one coral from this group that is overrepresented
    coral_to_move <- group_rows[
      ssids_assigned$sample_time[group_rows] == over_time
    ][1]
    this_tank <- ssids_assigned$tank[coral_to_move]
    # find a different coral in same tank that is currently in the missing timepoint
    swap_candidate <- which(
      ssids_assigned$tank == this_tank &
        ssids_assigned$sample_time == under_time &
        !(
          ssids_assigned$site == this_site &
            ssids_assigned$treatment == this_treatment
        )
    )
    # if no valid swap exists, stop trying this group
    if(length(swap_candidate) == 0) {
      warning(paste(
        "Could not fully fix:",
        this_site,
        this_treatment
      ))
      break
    }
    swap_candidate <- sample(swap_candidate, 1)
    # swap timepoints only
    ssids_assigned$sample_time[coral_to_move] <- under_time
    ssids_assigned$sample_time[swap_candidate] <- over_time
  }
}


## ------------------ CHECKS ------------------------------

table(ssids_assigned$sample_time, ssids_assigned$tank)

table(ssids_assigned$sample_time, ssids_assigned$site)

table(ssids_assigned$sample_time, ssids_assigned$treatment)

table(ssids_assigned$site, ssids_assigned$treatment)


results <- ssids_assigned %>%
  
  count(site, treatment, sample_time) %>%
  
  pivot_wider(
    
    names_from = sample_time,
    
    values_from = n,
    
    values_fill = 0
    
  )

results


# -------------------------- Clean table and export ---------------------------


final <- frags %>%
  select(-site_group, -priority, -id)

write.csv(final, "fragment_placement.csv", row.names = FALSE)

