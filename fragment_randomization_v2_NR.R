
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






# ------------------ Assigning time point sampling with swap optimizer ------------------



set.seed(123)

n_timepoints <- 5
timepoints <- paste0("T", 1:n_timepoints)
per_timepoint <- 5


# Function to front-load target counts
# n = 4  -> 1,1,1,1,0
# n = 5  -> 1,1,1,1,1
# n = 12 -> 3,3,2,2,2
make_time_targets <- function(n, n_timepoints = 5) {
  base <- floor(n / n_timepoints)
  extra <- n %% n_timepoints
  
  out <- rep(base, n_timepoints)
  
  if (extra > 0) {
    out[1:extra] <- out[1:extra] + 1
  }
  
  out
}

# Add row ID and group ID
ssids_work <- ssids_assigned %>%
  mutate(
    row_id = row_number(),
    site_treatment = paste(site, treatment, sep = "__")
  )

# ------------------ 1. Create ideal targets for site x treatment ------------------

target_df <- ssids_work %>%
  count(site_treatment, site, treatment, name = "n_total") %>%
  rowwise() %>%
  mutate(targets = list(make_time_targets(n_total, n_timepoints))) %>%
  ungroup() %>%
  unnest_wider(targets, names_sep = "_") %>%
  rename(
    T1 = targets_1,
    T2 = targets_2,
    T3 = targets_3,
    T4 = targets_4,
    T5 = targets_5
  )

target_long <- target_df %>%
  select(site_treatment, all_of(timepoints)) %>%
  pivot_longer(
    cols = all_of(timepoints),
    names_to = "sample_time",
    values_to = "target"
  )

# ------------------ 2. Initial assignment: exactly 5 per tank per timepoint ------------------

ssids_sampled <- ssids_work %>%
  group_by(tank) %>%
  mutate(
    sample_time = sample(rep(timepoints, each = per_timepoint))
  ) %>%
  ungroup()

# ------------------ 3. Objective function ------------------

calc_score <- function(dat, target_long) {
  
  current <- dat %>%
    count(site_treatment, sample_time, name = "current") %>%
    right_join(target_long, by = c("site_treatment", "sample_time")) %>%
    mutate(current = replace_na(current, 0))
  
  # Basic mismatch from ideal target
  base_score <- sum(abs(current$current - current$target))
  
  # Extra penalty when a group has 0 in an earlier timepoint
  # but has individuals in a later timepoint.
  zero_penalty <- current %>%
    mutate(time_num = as.integer(gsub("T", "", sample_time))) %>%
    group_by(site_treatment) %>%
    arrange(time_num, .by_group = TRUE) %>%
    summarise(
      penalty = sum(
        sapply(seq_len(n()), function(i) {
          if (current[i] == 0 && any(current[(i + 1):n()] > 0, na.rm = TRUE)) {
            10
          } else {
            0
          }
        })
      ),
      .groups = "drop"
    ) %>%
    summarise(total = sum(penalty)) %>%
    pull(total)
  
  base_score + zero_penalty
}

score <- calc_score(ssids_sampled, target_long)

score


# ------------------ 4. Swap optimizer ------------------

set.seed(123)

n_iter <- 220000 

score_history <- tibble(
  iteration = integer(),
  score = numeric(),
  best_score = numeric()
)

best_dat <- ssids_sampled
best_score <- score

for (i in seq_len(n_iter)) {
  
  # Temperature decreases over time
  temp <- 5 * (1 - i / n_iter)
  
  # Print progress
  if (i %% 1000 == 0) {
    message(
      "Iteration: ", i, " / ", n_iter,
      " | score = ", score,
      " | best score = ", best_score
    )
    
    score_history <- bind_rows(
      score_history,
      tibble(
        iteration = i,
        score = score,
        best_score = best_score
      )
    )
  }
  
  # Pick a random tank
  this_tank <- sample(unique(ssids_sampled$tank), 1)
  
  tank_rows <- ssids_sampled %>%
    filter(tank == this_tank) %>%
    pull(row_id)
  
  # Pick two corals in that tank
  swap_rows <- sample(tank_rows, 2)
  
  time_1 <- ssids_sampled$sample_time[ssids_sampled$row_id == swap_rows[1]]
  time_2 <- ssids_sampled$sample_time[ssids_sampled$row_id == swap_rows[2]]
  
  if (time_1 == time_2) next
  
  # Try swap
  dat_try <- ssids_sampled
  
  dat_try$sample_time[dat_try$row_id == swap_rows[1]] <- time_2
  dat_try$sample_time[dat_try$row_id == swap_rows[2]] <- time_1
  
  new_score <- calc_score(dat_try, target_long)
  
  delta <- new_score - score
  
  # Accept if better OR sometimes if worse
  accept <- delta <= 0 || runif(1) < exp(-delta / max(temp, 0.001))
  
  if (accept) {
    ssids_sampled <- dat_try
    score <- new_score
  }
  
  # Save best version ever found
  if (score < best_score) {
    best_score <- score
    best_dat <- ssids_sampled
  }
  
  if (best_score == 0) {
    message("Perfect score reached at iteration ", i)
    break
  }
}

# Use the best version found, not necessarily the final version
ssids_sampled <- best_dat
score <- best_score

score




ssids_assigned_sampled <- ssids_sampled %>%
  select(-row_id, -site_treatment)





# ------------------ Checks ------------------

# 1. Tank x timepoint should be exactly 5 everywhere
ssids_assigned_sampled %>%
  count(tank, sample_time) %>%
  pivot_wider(
    names_from = sample_time,
    values_from = n,
    values_fill = 0
  ) %>%
  arrange(tank)

# 2. Each tank should still have 25
ssids_assigned_sampled %>%
  count(tank)

# 3. Site x treatment x timepoint balance
results = ssids_assigned_sampled %>%
  count(site, treatment, sample_time) %>%
  pivot_wider(
    names_from = sample_time,
    values_from = n,
    values_fill = 0
  ) %>%
  arrange(site, treatment)

# 4. Dodge Island
ssids_assigned_sampled %>%
  filter(site == "Dodge Island") %>%
  count(site, treatment, sample_time) %>%
  pivot_wider(
    names_from = sample_time,
    values_from = n,
    values_fill = 0
  )

# 5. Hollywood
ssids_assigned_sampled %>%
  filter(site == "Hollywood") %>%
  count(site, treatment, sample_time) %>%
  pivot_wider(
    names_from = sample_time,
    values_from = n,
    values_fill = 0
  )


write_csv(ssids_assigned_sampled, "final_ssid_assignments.csv")
