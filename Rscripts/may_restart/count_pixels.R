# Suppose your raster is
r <- result_raster$SCANFI_1km

# 1. Count number of pixels per class
pixel_counts <- terra::freq(r)  # returns a data.frame with 'value' and 'count'

# 2. Compute ratios relative to the class with the most pixels
max_count <- max(pixel_counts$count)
pixel_counts$ratio <- pixel_counts$count / max_count

pixel_counts

# convert to numeric
confusion[[1]] <- confusion[[1]] %>%
  mutate(actual = as.numeric(actual),
         predicted = as.numeric(predicted))

# Summarize counts
actual_counts <- confusion[[1]] %>%
  count(actual, name = "count_actual")

pred_counts <- confusion[[1]] %>%
  count(predicted, name = "count_pred")

training <- tibble(class = as.numeric(y_int[-holdout_idx])) |> 
  group_by(class) |> 
  count() |>
  rename(count_train = n)

# Join with raster pixel counts
comparison <- pixel_counts %>%
  rename(class = value, pixel_count = count) %>%
  left_join(actual_counts, by = c("class" = "actual")) %>%
  left_join(pred_counts, by = c("class" = "predicted")) %>%
  left_join(training)

# Optional: ratios relative to max
comparison <- comparison %>%
  mutate(pixel_ratio = pixel_count / max(pixel_count, na.rm = TRUE),
         actual_ratio = count_actual / max(count_actual, na.rm = TRUE),
         pred_ratio = count_pred / max(count_pred, na.rm = TRUE),
         train_ratio = count_train / max(count_train, na.rm=TRUE))

comparison


  
