# ------------------------------------------------------------------------------
# Test Behaviour of recordLinkage()
# ------------------------------------------------------------------------------

devtools::load_all()

# Result recordLinkage is Row Order Dependent ----------------------------------
# depending on how the rows are sorted, the results differ

# test data (standard row order)
df <- data.frame(age = c(20, 20, 20, 25, 25),
                 gender = c("f", "f", "f", "m", "m"),
                 sens = c(1,1,2,3,4))

df_ano <- df
df_ano[4, "gender"] <- "f" # PRAM

df; df_ano

# record linkage
rL1 <- recordLinkage(df, df_ano,
                    key = c("age", "gender"),
                    strategy = "nearest",
                    return_matches = TRUE,
                    matching = "bijective")
rL1
rL1$per_record

# test data (row order changes in anonymized df)
df$id <- 1:5
df_ano$id <- 1:5
df_ano <- df_ano[c(3,2,1,4,5),] # change order of the first 3 rows


df; df_ano

rL2 <- recordLinkage(df, df_ano,
                       key = c("age", "gender"),
                       strategy = "nearest",
                       return_matches = TRUE,
                       truth = "id",
                       id = "id",
                       matching = "bijective")
rL2
rL2$per_record

# Edge Case with unrelated Data Frames -----------------------------------------
# extreme case that shows why row order matters

# test data
df <- data.frame(occupation = c("priest", "bartender", "pilot", "teacher", "scientist"),
                 degree = c("bachelor", "bachelor", "bachelor", "master", "master"))

df_ano <- data.frame(occupation = c("hunter", "sailor", "carpenter", "president", "waiter"),
                     degree = c("phd", "phd", "phd", "phd", "phd"))

df; df_ano

# Risk should be very low; any observed match would be due to chance alone.
rL <- recordLinkage(df, df_ano,
                    key = c("occupation", "degree"),
                    strategy = "nearest",
                    return_matches = TRUE,
                    matching = "bijective")

rL
rL$per_record

# Potential Solution -----------------------------------------------------------

# Random permutation of the rows of the anonymized dataset
# Assess the average risk for each permutation
# Show distribution of the different risk scores
